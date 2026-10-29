from __future__ import annotations

import csv
import json
from collections import Counter
from collections.abc import Iterable, Mapping
from pathlib import Path

from reaxnet.models import (
    AnalysisResult,
    BondKey,
    FrameState,
    MoleculeRecord,
    ReactionEvent,
)


SPECIES_FIELDS = [
    "structure_id",
    "timestep",
    "time_ps",
    "species_key",
    "formula",
    "smiles",
    "inchi",
    "inchikey",
    "count",
]
EVENT_FIELDS = [
    "event_id",
    "previous_timestep",
    "timestep",
    "time_ps",
    "event_type",
    "reaction",
    "reactants",
    "products",
    "reactant_keys",
    "product_keys",
    "broken_bonds",
    "formed_bonds",
    "changed_atom_ids",
    "involved_atom_ids",
    "accepted",
    "filter_reason",
    "filter_group",
]
REGISTRY_FIELDS = [
    "structure_id",
    "species_key",
    "formula",
    "atom_count",
    "bond_count",
    "example_atom_ids",
    "bond_signature",
    "smiles",
    "inchi",
    "inchikey",
    "status",
    "message",
]


def _number(value: float) -> str:
    return format(value, ".12g")


def _join(values: Iterable[object]) -> str:
    return ";".join(str(value) for value in values)


def _bond_list(bonds: tuple[BondKey, ...]) -> str:
    return ";".join(f"{atom1}-{atom2}" for atom1, atom2 in bonds)


def _label(key: str, registry: Mapping[str, MoleculeRecord]) -> str:
    molecule = registry[key]
    return molecule.smiles or molecule.formula


def _reaction_side(keys: tuple[str, ...], registry: Mapping[str, MoleculeRecord]) -> str:
    return " + ".join(_label(key, registry) for key in keys)


def _write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def species_rows_for_state(
    state: FrameState,
    registry: Mapping[str, MoleculeRecord],
    timestep_scale: float,
) -> Iterable[dict[str, object]]:
    counts = Counter(molecule.species_key for molecule in state.molecules)
    for key, count in sorted(counts.items()):
        molecule = registry[key]
        yield {
            "structure_id": molecule.structure_id,
            "timestep": state.timestep,
            "time_ps": _number(state.timestep * timestep_scale),
            "species_key": key,
            "formula": molecule.formula,
            "smiles": molecule.smiles,
            "inchi": molecule.inchi,
            "inchikey": molecule.inchikey,
            "count": count,
        }


def event_to_row(
    event: ReactionEvent,
    registry: Mapping[str, MoleculeRecord],
) -> dict[str, object]:
    reactants = _reaction_side(event.reactant_keys, registry)
    products = _reaction_side(event.product_keys, registry)
    return {
        "event_id": event.event_id,
        "previous_timestep": event.previous_timestep,
        "timestep": event.timestep,
        "time_ps": _number(event.time_ps),
        "event_type": event.event_type,
        "reaction": f"{reactants} -> {products}",
        "reactants": reactants,
        "products": products,
        "reactant_keys": _join(event.reactant_keys),
        "product_keys": _join(event.product_keys),
        "broken_bonds": _bond_list(event.broken_bonds),
        "formed_bonds": _bond_list(event.formed_bonds),
        "changed_atom_ids": _join(event.changed_atom_ids),
        "involved_atom_ids": _join(event.involved_atom_ids),
        "accepted": str(event.accepted).lower(),
        "filter_reason": event.filter_reason,
        "filter_group": event.filter_group,
    }


def event_to_payload(event: ReactionEvent, registry: Mapping[str, MoleculeRecord]) -> dict:
    """Finalize mutable filter/ID fields while keeping original frame-local evidence."""
    row = event_to_row(event, registry)
    def structure_side(molecules):
        counts = Counter(m.get("structure_id") or m["species_key"] for m in molecules)
        return " + ".join((f"{count} " if count != 1 else "") + key
                          for key,count in sorted(counts.items()))
    reactants = event.evidence.get("reactants", [])
    products = event.evidence.get("products", [])
    reaction = (structure_side(reactants) + " -> " + structure_side(products)
                if reactants and products else row["reaction"])
    return {
        **event.evidence,
        "event_id": event.event_id,
        "previous_timestep": event.previous_timestep, "timestep": event.timestep,
        "time_start_ps": event.evidence.get("time_start_ps"),
        "time_end_ps": event.evidence.get("time_end_ps", event.time_ps),
        "reaction": reaction, "reaction_identifiers": row["reaction"],
        "event_type": event.event_type,
        "broken_bonds": [list(edge) for edge in event.broken_bonds],
        "formed_bonds": [list(edge) for edge in event.formed_bonds],
        "changed_atom_ids": list(event.changed_atom_ids),
        "involved_atom_ids": list(event.involved_atom_ids),
        "accepted": event.accepted, "filter_reason": event.filter_reason,
        "filter_group": event.filter_group,
    }


def registry_to_row(key: str, molecule: MoleculeRecord) -> dict[str, object]:
    return {
        "structure_id": molecule.structure_id,
        "species_key": key,
        "formula": molecule.formula,
        "atom_count": molecule.atom_count,
        "bond_count": molecule.bond_count,
        "example_atom_ids": _join(molecule.atom_ids),
        "bond_signature": molecule.bond_signature,
        "smiles": molecule.smiles,
        "inchi": molecule.inchi,
        "inchikey": molecule.inchikey,
        "status": molecule.status,
        "message": molecule.message,
    }


def write_outputs(result: AnalysisResult, output_dir: str | Path) -> dict[str, Path]:
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    paths = {
        "species_timeseries": target / "species_timeseries.csv",
        "reaction_events": target / "reaction_events.csv",
        "molecule_registry": target / "molecule_registry.csv",
        "reaction_network": target / "reaction_network.json",
    }

    _write_csv(
        paths["species_timeseries"],
        SPECIES_FIELDS,
        (
            row
            for state in result.states
            for row in species_rows_for_state(
                state,
                result.registry,
                result.config.timestep_scale,
            )
        ),
    )

    _write_csv(
        paths["reaction_events"],
        EVENT_FIELDS,
        (event_to_row(event, result.registry) for event in result.events),
    )

    _write_csv(
        paths["molecule_registry"],
        REGISTRY_FIELDS,
        (registry_to_row(key, molecule) for key, molecule in sorted(result.registry.items())),
    )

    with paths["reaction_network"].open("w", encoding="utf-8") as handle:
        json.dump(result.network, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    return paths
