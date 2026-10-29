"""Bounded-memory ReaxFF analysis for multi-gigabyte bonds trajectories."""

from __future__ import annotations

import csv
import json
from collections import Counter
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path
from time import monotonic
from typing import Any

from reaxnet.analysis import (
    AnalyzerConfig,
    _events_between,
    _identity_rank,
    filter_multistep_recrossing,
)
from reaxnet.chemistry import build_frame_state, TopologyRegistry
from reaxnet.connectivity import ConnectivityConfig, ConnectivityTracker
from reaxnet.counts import CountWriter
from reaxnet.event_index import EventIndex
from reaxnet.identity import IdentityResult
from reaxnet.models import FrameState, MoleculeRecord, ReactionEvent
from reaxnet.outputs import (
    EVENT_FIELDS,
    REGISTRY_FIELDS,
    SPECIES_FIELDS,
    event_to_row,
    registry_to_row,
    species_rows_for_state,
)
from reaxnet.parser import parse_bonds_lines


ProgressCallback = Callable[[dict[str, int | float]], None]


@dataclass(frozen=True)
class StreamingAnalysisResult:
    paths: dict[str, Path]
    network: dict[str, Any]
    registry: dict[str, MoleculeRecord]
    summary: dict[str, int | bool | None]


def analyze_file_streaming(
    bonds_path: str | Path,
    element_mapping: Mapping[int, str],
    output_dir: str | Path,
    config: AnalyzerConfig | None = None,
    *,
    progress_callback: ProgressCallback | None = None,
) -> StreamingAnalysisResult:
    """Analyze a trajectory while retaining only the active recrossing window."""

    config = config or AnalyzerConfig()
    source = Path(bonds_path).resolve()
    source_size = source.stat().st_size
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    paths = {
        "species_timeseries": target / "species_timeseries.csv",
        "reaction_events": target / "reaction_events.csv",
        "molecule_registry": target / "molecule_registry.csv",
        "reaction_network": target / "reaction_network.json",
        "event_index": target / "event_index.sqlite",
    }
    count_writer = CountWriter(
        target, element_mapping, config.timestep_scale,
        source_path=source, expected_interval=config.expected_interval,
    )
    paths.update(count_writer.paths)

    registry: dict[str, MoleculeRecord] = {}
    identity_cache: dict[object, IdentityResult] = {}
    tracker = ConnectivityTracker(element_mapping, ConnectivityConfig(expected_interval=config.expected_interval))
    topologies = TopologyRegistry()
    frame_occurrences: Counter[str] = Counter()
    total_observations: Counter[str] = Counter()
    max_counts: Counter[str] = Counter()
    aggregates: dict[tuple[tuple[str, ...], tuple[str, ...]], dict[str, Any]] = {}
    state_buffer: list[FrameState] = []
    pending_events: list[ReactionEvent] = []
    previous_state: FrameState | None = None
    previous_timestep: int | None = None
    first_timestep: int | None = None
    last_timestep: int | None = None
    frame_count = 0
    event_number = 0
    recrossing_group_number = 0
    detected_event_count = 0
    accepted_event_count = 0
    filtered_event_count = 0
    transition_bond_count = 0
    peak_buffered_states = 0
    peak_pending_events = 0
    last_progress_at = 0.0

    with (
        count_writer,
        EventIndex(paths["event_index"]) as event_index,
        source.open("r", encoding="utf-8", errors="replace") as source_handle,
        paths["species_timeseries"].open("w", encoding="utf-8", newline="") as species_handle,
        paths["reaction_events"].open("w", encoding="utf-8", newline="") as event_handle,
    ):
        species_writer = csv.DictWriter(species_handle, fieldnames=SPECIES_FIELDS)
        event_writer = csv.DictWriter(event_handle, fieldnames=EVENT_FIELDS)
        species_writer.writeheader()
        event_writer.writeheader()
        lines = iter(source_handle.readline, "")

        for frame in parse_bonds_lines(lines):
            if previous_timestep is not None and frame.timestep <= previous_timestep:
                raise ValueError(
                    f"Timesteps must be strictly increasing: {previous_timestep} then {frame.timestep}"
                )
            previous_timestep = frame.timestep
            first_timestep = frame.timestep if first_timestep is None else first_timestep
            last_timestep = frame.timestep
            observed = tracker.process(frame)
            count_writer.write(frame, observed)
            state = build_frame_state(
                frame,
                element_mapping,
                generate_identifiers=config.generate_identifiers,
                identity_cache=identity_cache,
                connectivity_frame=observed,
                topology_registry=topologies,
            )
            frame_count += 1
            transition_bond_count += len(state.transition_bonds)
            state_buffer.append(state)

            for molecule in state.molecules:
                existing = registry.get(molecule.species_key)
                if existing is None or _identity_rank(molecule.status) > _identity_rank(existing.status):
                    registry[molecule.species_key] = molecule

            counts = Counter(molecule.species_key for molecule in state.molecules)
            for key, count in counts.items():
                frame_occurrences[key] += 1
                total_observations[key] += count
                max_counts[key] = max(max_counts[key], count)
            species_writer.writerows(
                species_rows_for_state(state, registry, config.timestep_scale)
            )

            if previous_state is not None:
                new_events = _events_between(previous_state, state, config.timestep_scale)
                new_events.sort(key=lambda event: event.changed_atom_ids)
                for event in new_events:
                    event_number += 1
                    event.event_id = f"E{event_number:06d}"
                pending_events.extend(new_events)

            if config.recrossing_window > 0 and len(pending_events) > 1:
                recrossing_group_number = filter_multistep_recrossing(
                    state_buffer,
                    pending_events,
                    config.recrossing_window,
                    recrossing_group_number,
                )

            cutoff_time = state.timestep * config.timestep_scale - config.recrossing_window
            if config.recrossing_window <= 0:
                ready = pending_events
                pending_events = []
            else:
                ready = [
                    event
                    for event in pending_events
                    if event.time_ps <= cutoff_time + 1e-12
                ]
                pending_events = [
                    event
                    for event in pending_events
                    if event.time_ps > cutoff_time + 1e-12
                ]

            for event in ready:
                event_writer.writerow(event_to_row(event, registry))
                event_index.write(event, registry)
                detected_event_count += 1
                if event.accepted:
                    accepted_event_count += 1
                    _aggregate_event(aggregates, event)
                else:
                    filtered_event_count += 1

            if pending_events:
                earliest_required = min(event.previous_timestep for event in pending_events)
                state_buffer = [
                    buffered_state
                    for buffered_state in state_buffer
                    if buffered_state.timestep >= earliest_required
                ]
            else:
                state_buffer = [state]

            peak_buffered_states = max(peak_buffered_states, len(state_buffer))
            peak_pending_events = max(peak_pending_events, len(pending_events))
            previous_state = state

            now = monotonic()
            if progress_callback and (frame_count == 1 or now - last_progress_at >= 0.5):
                progress_callback(
                    {
                        "frames_processed": frame_count,
                        "bytes_processed": min(_safe_tell(source_handle), source_size),
                        "bytes_total": source_size,
                    }
                )
                last_progress_at = now

        if frame_count == 0:
            raise ValueError("The bonds trajectory contains no timestep frames")

        for event in pending_events:
            event_writer.writerow(event_to_row(event, registry))
            event_index.write(event, registry)
            detected_event_count += 1
            if event.accepted:
                accepted_event_count += 1
                _aggregate_event(aggregates, event)
            else:
                filtered_event_count += 1

    network = _build_network(
        registry=registry,
        aggregates=aggregates,
        frame_occurrences=frame_occurrences,
        total_observations=total_observations,
        max_counts=max_counts,
        config=config,
        frame_count=frame_count,
        first_timestep=first_timestep,
        last_timestep=last_timestep,
        detected_event_count=detected_event_count,
        accepted_event_count=accepted_event_count,
        filtered_event_count=filtered_event_count,
        source_size=source_size,
        peak_buffered_states=peak_buffered_states,
        peak_pending_events=peak_pending_events,
    )
    network["metadata"]["observed_interval"] = tracker.interval
    with paths["molecule_registry"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REGISTRY_FIELDS)
        writer.writeheader()
        writer.writerows(
            registry_to_row(key, molecule)
            for key, molecule in sorted(registry.items())
        )
    with paths["reaction_network"].open("w", encoding="utf-8") as handle:
        json.dump(network, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    if progress_callback:
        progress_callback(
            {
                "frames_processed": frame_count,
                "bytes_processed": source_size,
                "bytes_total": source_size,
            }
        )

    summary: dict[str, int | bool | None] = {
        "frame_count": frame_count,
        "species_count": len(registry),
        "detected_event_count": detected_event_count,
        "accepted_event_count": accepted_event_count,
        "filtered_event_count": filtered_event_count,
        "reaction_count": len(network["reactions"]),
        "transition_bond_count": transition_bond_count,
        "source_size_bytes": source_size,
        "peak_buffered_states": peak_buffered_states,
        "peak_pending_events": peak_pending_events,
        "streaming_analysis": True,
        "observed_interval": tracker.interval,
    }
    return StreamingAnalysisResult(
        paths=paths,
        network=network,
        registry=registry,
        summary=summary,
    )


def _safe_tell(handle) -> int:
    try:
        return int(handle.tell())
    except (OSError, ValueError):
        return 0


def _aggregate_event(
    aggregates: dict[tuple[tuple[str, ...], tuple[str, ...]], dict[str, Any]],
    event: ReactionEvent,
) -> None:
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


def _side_payload(
    keys: tuple[str, ...],
    registry: Mapping[str, MoleculeRecord],
) -> list[dict[str, Any]]:
    return [
        {
            "species_key": key,
            "formula": registry[key].formula,
            "smiles": registry[key].smiles,
            "stoichiometry": count,
        }
        for key, count in sorted(Counter(keys).items())
    ]


def _build_network(
    *,
    registry: Mapping[str, MoleculeRecord],
    aggregates: Mapping[tuple[tuple[str, ...], tuple[str, ...]], dict[str, Any]],
    frame_occurrences: Counter[str],
    total_observations: Counter[str],
    max_counts: Counter[str],
    config: AnalyzerConfig,
    frame_count: int,
    first_timestep: int | None,
    last_timestep: int | None,
    detected_event_count: int,
    accepted_event_count: int,
    filtered_event_count: int,
    source_size: int,
    peak_buffered_states: int,
    peak_pending_events: int,
) -> dict[str, Any]:
    species = [
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
        for key, molecule in sorted(registry.items())
    ]
    reactions = []
    sorted_aggregates = sorted(
        aggregates.values(),
        key=lambda item: (item["first_timestep"], item["reactant_keys"]),
    )
    for index, aggregate in enumerate(sorted_aggregates, start=1):
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
    return {
        "metadata": {
            "frame_count": frame_count,
            "first_timestep": first_timestep,
            "last_timestep": last_timestep,
            "connectivity_policy": "file_neighbor_list",
            "timestep_scale_ps": config.timestep_scale,
            "expected_interval": config.expected_interval,
            "recrossing_window_ps": config.recrossing_window,
            "detected_event_count": detected_event_count,
            "accepted_event_count": accepted_event_count,
            "filtered_event_count": filtered_event_count,
            "streaming_analysis": True,
            "source_size_bytes": source_size,
            "peak_buffered_states": peak_buffered_states,
            "peak_pending_events": peak_pending_events,
        },
        "species": species,
        "reactions": reactions,
    }
