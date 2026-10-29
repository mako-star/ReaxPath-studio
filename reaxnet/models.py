from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


BondKey = tuple[int, int]


@dataclass(frozen=True)
class AtomRecord:
    atom_id: int
    atom_type: int
    neighbors: tuple[int, ...]
    molecule_id: int
    bond_orders: tuple[float, ...]
    abo: float
    lone_pairs: float
    charge: float


@dataclass(frozen=True)
class Frame:
    timestep: int
    atoms: dict[int, AtomRecord]
    expected_atom_count: int | None = None
    output_cutoff: float | None = None
    source_line_start: int | None = None
    source_line_end: int | None = None


@dataclass(frozen=True)
class Bond:
    atom1: int
    atom2: int
    raw_order: float
    effective_order: float

    @property
    def key(self) -> BondKey:
        return (self.atom1, self.atom2)


@dataclass(frozen=True)
class MoleculeRecord:
    atom_ids: tuple[int, ...]
    formula: str
    species_key: str
    bond_signature: str
    atom_count: int
    bond_count: int
    smiles: str = ""
    inchi: str = ""
    inchikey: str = ""
    status: str = "unavailable"
    message: str = ""
    instance_id: str = ""
    structure_id: str = ""
    atom_elements: tuple[str, ...] = ()
    connections: tuple[BondKey, ...] = ()


@dataclass(frozen=True)
class FrameState:
    timestep: int
    bonds: dict[BondKey, Bond]
    transition_bonds: tuple[Bond, ...]
    molecules: tuple[MoleculeRecord, ...]
    atom_to_molecule: dict[int, int]
    symbols: dict[int, str] = field(default_factory=dict)


@dataclass
class ReactionEvent:
    event_id: str
    previous_timestep: int
    timestep: int
    time_ps: float
    event_type: str
    reactant_keys: tuple[str, ...]
    product_keys: tuple[str, ...]
    reactant_formulas: tuple[str, ...]
    product_formulas: tuple[str, ...]
    broken_bonds: tuple[BondKey, ...]
    formed_bonds: tuple[BondKey, ...]
    changed_atom_ids: tuple[int, ...]
    involved_atom_ids: tuple[int, ...]
    accepted: bool = True
    filter_reason: str = ""
    filter_group: str = ""
    evidence: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class AnalysisResult:
    states: tuple[FrameState, ...]
    events: tuple[ReactionEvent, ...]
    registry: dict[str, MoleculeRecord]
    network: dict[str, Any]
    config: Any = field(repr=False)
