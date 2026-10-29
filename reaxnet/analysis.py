from __future__ import annotations

from collections import Counter, defaultdict, deque
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from reaxnet.chemistry import build_frame_state, TopologyRegistry
from reaxnet.connectivity import ConnectivityConfig, ConnectivityTracker
import math
from reaxnet.identity import IdentityResult
from reaxnet.models import (
    AnalysisResult,
    BondKey,
    Frame,
    FrameState,
    MoleculeRecord,
    ReactionEvent,
)
from reaxnet.parser import parse_bonds_file


@dataclass(frozen=True)
class AnalyzerConfig:
    timestep_scale: float = 1.0
    recrossing_window: float = 0.0
    generate_identifiers: bool = True
    expected_interval: int | None = None

    def __post_init__(self) -> None:
        ConnectivityConfig(expected_interval=self.expected_interval)
        if not math.isfinite(self.timestep_scale) or self.timestep_scale <= 0:
            raise ValueError("timestep_scale must be positive")
        if not math.isfinite(self.recrossing_window) or self.recrossing_window < 0:
            raise ValueError("recrossing_window cannot be negative")


def analyze_file(
    bonds_path: str | Path,
    element_mapping: Mapping[int, str],
    config: AnalyzerConfig | None = None,
) -> AnalysisResult:
    return analyze_frames(parse_bonds_file(bonds_path), element_mapping, config)


def analyze_frames(
    frames: Iterable[Frame],
    element_mapping: Mapping[int, str],
    config: AnalyzerConfig | None = None,
) -> AnalysisResult:
    config = config or AnalyzerConfig()
    states: list[FrameState] = []
    registry: dict[str, MoleculeRecord] = {}
    identity_cache: dict[object, IdentityResult] = {}
    previous_timestep: int | None = None
    tracker = ConnectivityTracker(element_mapping, ConnectivityConfig(expected_interval=config.expected_interval))
    topologies = TopologyRegistry()

    for frame in frames:
        if previous_timestep is not None and frame.timestep <= previous_timestep:
            raise ValueError(
                f"Timesteps must be strictly increasing: {previous_timestep} then {frame.timestep}"
            )
        previous_timestep = frame.timestep
        observed = tracker.process(frame)
        state = build_frame_state(
            frame,
            element_mapping,
            generate_identifiers=config.generate_identifiers,
            identity_cache=identity_cache,
            connectivity_frame=observed,
            topology_registry=topologies,
        )
        states.append(state)
        for molecule in state.molecules:
            existing = registry.get(molecule.species_key)
            if existing is None or _identity_rank(molecule.status) > _identity_rank(existing.status):
                registry[molecule.species_key] = molecule

    if not states:
        raise ValueError("The bonds trajectory contains no timestep frames")

    events: list[ReactionEvent] = []
    for old_state, new_state in zip(states, states[1:], strict=False):
        events.extend(_events_between(old_state, new_state, config.timestep_scale))
    events.sort(key=lambda event: (event.timestep, event.changed_atom_ids))
    for index, event in enumerate(events, start=1):
        event.event_id = f"E{index:06d}"

    if config.recrossing_window > 0:
        filter_multistep_recrossing(states, events, config.recrossing_window)

    network = build_reaction_network(states, events, registry, config)
    return AnalysisResult(
        states=tuple(states),
        events=tuple(events),
        registry=registry,
        network=network,
        config=config,
    )


def _identity_rank(status: str) -> int:
    return {"ok": 4, "fallback": 3, "partial": 2, "disabled": 1, "unavailable": 1}.get(
        status, 0
    )


def _events_between(old: FrameState, new: FrameState, timestep_scale: float) -> list[ReactionEvent]:
    if set(old.atom_to_molecule) != set(new.atom_to_molecule):
        raise ValueError(
            f"Atom ids changed between timesteps {old.timestep} and {new.timestep}; "
            "variable-atom trajectories are not supported"
        )

    broken = set(old.bonds) - set(new.bonds)
    formed = set(new.bonds) - set(old.bonds)
    if not broken and not formed:
        return []

    changed_atoms = {atom_id for bond in broken | formed for atom_id in bond}
    graph: dict[tuple[str, int], set[tuple[str, int]]] = defaultdict(set)
    for atom_id in old.atom_to_molecule:
        old_node = ("old", old.atom_to_molecule[atom_id])
        new_node = ("new", new.atom_to_molecule[atom_id])
        graph[old_node].add(new_node)
        graph[new_node].add(old_node)

    components: list[set[tuple[str, int]]] = []
    unseen = set(graph)
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        queue = deque([start])
        component = {start}
        while queue:
            node = queue.popleft()
            for neighbor in graph[node]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    component.add(neighbor)
                    queue.append(neighbor)
        components.append(component)

    events: list[ReactionEvent] = []
    for component in components:
        if not any(("old", old.atom_to_molecule[a]) in component for a in changed_atoms):
            continue
        old_indexes = sorted(index for side, index in component if side == "old")
        new_indexes = sorted(index for side, index in component if side == "new")
        reactants = sorted(
            (old.molecules[index] for index in old_indexes),
            key=lambda molecule: (molecule.species_key, molecule.atom_ids),
        )
        products = sorted(
            (new.molecules[index] for index in new_indexes),
            key=lambda molecule: (molecule.species_key, molecule.atom_ids),
        )
        component_changed_atoms = {
            atom_id
            for atom_id in changed_atoms
            if old.atom_to_molecule[atom_id] in old_indexes
            or new.atom_to_molecule[atom_id] in new_indexes
        }
        component_broken = tuple(
            sorted(bond for bond in broken if bond[0] in component_changed_atoms)
        )
        component_formed = tuple(
            sorted(bond for bond in formed if bond[0] in component_changed_atoms)
        )
        involved_atoms = sorted(
            {atom_id for molecule in reactants + products for atom_id in molecule.atom_ids}
        )
        old_ids = Counter(a for m in reactants for a in m.atom_ids)
        new_ids = Counter(a for m in products for a in m.atom_ids)
        if old_ids != new_ids or any(n != 1 for n in old_ids.values()):
            raise ValueError("Event does not conserve atom IDs exactly once on each side")
        old_elements = Counter(old.symbols[a] for a in old_ids)
        new_elements = Counter(new.symbols[a] for a in new_ids)
        if old_elements != new_elements:
            raise ValueError("Event does not conserve elements")
        evidence = {
            "time_start_ps": old.timestep * timestep_scale,
            "time_end_ps": new.timestep * timestep_scale,
            "reactants": [_molecule_evidence(m) for m in reactants],
            "products": [_molecule_evidence(m) for m in products],
            "conservation": {"atom_ids": True, "elements": True},
            "interpretation": "candidate_changes_between_saved_frames",
        }
        events.append(
            ReactionEvent(
                event_id="",
                previous_timestep=old.timestep,
                timestep=new.timestep,
                time_ps=new.timestep * timestep_scale,
                event_type=_classify_event(len(reactants), len(products)),
                reactant_keys=tuple(molecule.species_key for molecule in reactants),
                product_keys=tuple(molecule.species_key for molecule in products),
                reactant_formulas=tuple(molecule.formula for molecule in reactants),
                product_formulas=tuple(molecule.formula for molecule in products),
                broken_bonds=component_broken,
                formed_bonds=component_formed,
                changed_atom_ids=tuple(sorted(component_changed_atoms)),
                involved_atom_ids=tuple(involved_atoms),
                evidence=evidence,
            )
        )
    return events


def _molecule_evidence(molecule):
    return {
        "instance_id": molecule.instance_id,
        "species_key": molecule.species_key,
        "structure_id": molecule.structure_id,
        "formula": molecule.formula,
        "atoms": [{"atom_id": a, "element": e} for a,e in zip(molecule.atom_ids,molecule.atom_elements)],
        "connections": [list(edge) for edge in molecule.connections],
        "smiles": molecule.smiles, "inchi": molecule.inchi, "inchikey": molecule.inchikey,
        "identity_status": molecule.status, "message": molecule.message,
    }


def _classify_event(reactant_count: int, product_count: int) -> str:
    if reactant_count > 1 and product_count == 1:
        return "association"
    if reactant_count == 1 and product_count > 1:
        return "dissociation"
    if reactant_count == 1 and product_count == 1:
        return "isomerization"
    if reactant_count > 1 and product_count > 1:
        return "exchange"
    return "rearrangement"


def filter_multistep_recrossing(
    states: list[FrameState],
    events: list[ReactionEvent],
    window_ps: float,
    group_number: int = 0,
) -> int:
    """Filter atom-connected event chains that return to their initial bond state."""

    state_by_timestep = {state.timestep: state for state in states}
    for start_index, start in enumerate(events):
        if not start.accepted:
            continue
        initial_state = state_by_timestep[start.previous_timestep]

        for end_index in range(start_index + 1, len(events)):
            end = events[end_index]
            if end.time_ps - start.time_ps > window_ps + 1e-12:
                break
            if not end.accepted:
                continue

            included = _connected_event_closure(events, start_index, end_index)
            if end_index not in included or len(included) < 2:
                continue
            atom_ids = {
                atom_id
                for event_index in included
                for atom_id in events[event_index].involved_atom_ids
            }
            final_state = state_by_timestep[end.timestep]
            if _bond_state_signature(initial_state, atom_ids) != _bond_state_signature(
                final_state, atom_ids
            ):
                continue

            included_events = [events[index] for index in sorted(included)]
            exact_reverse = (
                len(included_events) == 2
                and included_events[0].reactant_keys == included_events[1].product_keys
                and included_events[0].product_keys == included_events[1].reactant_keys
            )
            reason = (
                "exact_reverse_within_window"
                if exact_reverse
                else "multistep_state_return_within_window"
            )
            group_number += 1
            group_id = f"RC{group_number:06d}"
            for event in included_events:
                event.accepted = False
                event.filter_reason = reason
                event.filter_group = group_id
            break
    return group_number


def _connected_event_closure(
    events: list[ReactionEvent],
    start_index: int,
    end_index: int,
) -> set[int]:
    included = {start_index}
    atom_ids = set(events[start_index].involved_atom_ids)
    changed = True
    while changed:
        changed = False
        for event_index in range(start_index + 1, end_index + 1):
            event = events[event_index]
            if not event.accepted or event_index in included:
                continue
            event_atoms = set(event.involved_atom_ids)
            if atom_ids & event_atoms:
                included.add(event_index)
                atom_ids.update(event_atoms)
                changed = True
    return included


def _bond_state_signature(
    state: FrameState,
    atom_ids: set[int],
) -> tuple[tuple[int, ...], tuple[BondKey, ...]]:
    selected_atoms = tuple(sorted(atom_ids))
    selected_bonds = tuple(
        sorted(
            bond
            for bond in state.bonds
            if bond[0] in atom_ids and bond[1] in atom_ids
        )
    )
    return selected_atoms, selected_bonds


def _side_payload(keys: tuple[str, ...], registry: Mapping[str, MoleculeRecord]) -> list[dict[str, Any]]:
    return [
        {
            "species_key": key,
            "formula": registry[key].formula,
            "smiles": registry[key].smiles,
            "stoichiometry": count,
        }
        for key, count in sorted(Counter(keys).items())
    ]


def build_reaction_network(
    states: list[FrameState],
    events: list[ReactionEvent],
    registry: Mapping[str, MoleculeRecord],
    config: AnalyzerConfig,
) -> dict[str, Any]:
    frame_occurrences = Counter()
    total_observations = Counter()
    max_counts = Counter()
    for state in states:
        counts = Counter(molecule.species_key for molecule in state.molecules)
        for key, count in counts.items():
            frame_occurrences[key] += 1
            total_observations[key] += count
            max_counts[key] = max(max_counts[key], count)

    species = []
    for key, molecule in sorted(registry.items()):
        species.append(
            {
                "id": key,
                "label": molecule.smiles or molecule.formula,
                "formula": molecule.formula,
                "smiles": molecule.smiles,
                "inchi": molecule.inchi,
                "inchikey": molecule.inchikey,
                "status": molecule.status,
                "frame_occurrences": frame_occurrences[key],
                "total_observations": total_observations[key],
                "max_count": max_counts[key],
            }
        )

    aggregates: dict[tuple[tuple[str, ...], tuple[str, ...]], dict[str, Any]] = {}
    for event in (event for event in events if event.accepted):
        side1, side2 = event.reactant_keys, event.product_keys
        aggregate_key = (side1, side2) if side1 <= side2 else (side2, side1)
        aggregate = aggregates.get(aggregate_key)
        if aggregate is None:
            aggregate = {
                "reactant_keys": side1,
                "product_keys": side2,
                "forward_count": 0,
                "reverse_count": 0,
                "first_timestep": event.timestep,
                "last_timestep": event.timestep,
                "first_time_ps": event.time_ps,
                "last_time_ps": event.time_ps,
                "event_types": set(),
            }
            aggregates[aggregate_key] = aggregate
        if side1 == aggregate["reactant_keys"] and side2 == aggregate["product_keys"]:
            aggregate["forward_count"] += 1
        else:
            aggregate["reverse_count"] += 1
        aggregate["last_timestep"] = event.timestep
        aggregate["last_time_ps"] = event.time_ps
        aggregate["event_types"].add(event.event_type)

    reactions = []
    for index, aggregate in enumerate(
        sorted(aggregates.values(), key=lambda item: (item["first_timestep"], item["reactant_keys"])),
        start=1,
    ):
        forward = aggregate["forward_count"]
        reverse = aggregate["reverse_count"]
        reactions.append(
            {
                "id": f"R{index:06d}",
                "reactants": _side_payload(aggregate["reactant_keys"], registry),
                "products": _side_payload(aggregate["product_keys"], registry),
                "event_types": sorted(aggregate["event_types"]),
                "forward_count": forward,
                "reverse_count": reverse,
                "event_count": forward + reverse,
                "net_count": forward - reverse,
                "first_timestep": aggregate["first_timestep"],
                "last_timestep": aggregate["last_timestep"],
                "first_time_ps": aggregate["first_time_ps"],
                "last_time_ps": aggregate["last_time_ps"],
            }
        )

    accepted_count = sum(event.accepted for event in events)
    return {
        "metadata": {
            "frame_count": len(states),
            "first_timestep": states[0].timestep,
            "last_timestep": states[-1].timestep,
            "connectivity_policy": "file_neighbor_list",
            "timestep_scale_ps": config.timestep_scale,
            "expected_interval": config.expected_interval,
            "recrossing_window_ps": config.recrossing_window,
            "detected_event_count": len(events),
            "accepted_event_count": accepted_count,
            "filtered_event_count": len(events) - accepted_count,
        },
        "species": species,
        "reactions": reactions,
    }
