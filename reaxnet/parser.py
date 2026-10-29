from __future__ import annotations

import re
import math
from collections.abc import Iterable, Iterator
from pathlib import Path

from reaxnet.models import AtomRecord, Frame


TIMESTEP_RE = re.compile(r"^#\s*Timestep\s+(-?\d+)\s*$", re.IGNORECASE)
ATOM_COUNT_RE = re.compile(r"^#\s*Number\s+of\s+particles\s+(\d+)\s*$", re.IGNORECASE)

OUTPUT_CUTOFF_RE = re.compile(r"coarse bond order cutoff\s+(\S+)", re.IGNORECASE)

PERIODIC_SYMBOLS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
}


class BondsParseError(ValueError):
    """Raised when a ReaxFF bonds trajectory is structurally invalid."""


def canonical_symbol(value: str) -> str:
    symbol = value.strip().lower().capitalize()
    if symbol not in PERIODIC_SYMBOLS:
        raise ValueError(f"Unknown element symbol: {value!r}")
    return symbol


def parse_element_mapping(value: str) -> dict[int, str]:
    """Parse a LAMMPS atom-type mapping such as ``1:C,2:H,3:O``."""

    mapping: dict[int, str] = {}
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            raw_type, raw_symbol = item.split(":", 1)
            atom_type = int(raw_type.strip())
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid element mapping item {item!r}; expected TYPE:ELEMENT"
            ) from exc
        if atom_type <= 0:
            raise ValueError(f"LAMMPS atom type must be positive: {atom_type}")
        if atom_type in mapping:
            raise ValueError(f"LAMMPS atom type {atom_type} is mapped more than once")
        mapping[atom_type] = canonical_symbol(raw_symbol)
    if not mapping:
        raise ValueError("Element mapping is empty")
    return mapping


def parse_symbol_values(value: str, *, integer: bool = False) -> dict[str, float | int]:
    """Parse element-specific overrides such as ``N:3,S:6``."""

    parsed: dict[str, float | int] = {}
    if not value.strip():
        return parsed
    for item in value.split(","):
        try:
            raw_symbol, raw_number = item.split(":", 1)
            symbol = canonical_symbol(raw_symbol)
            number = int(raw_number) if integer else float(raw_number)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"Invalid override item {item!r}; expected ELEMENT:NUMBER"
            ) from exc
        if number < 0:
            raise ValueError(f"Override for {symbol} cannot be negative")
        parsed[symbol] = number
    return parsed


def parse_bonds_file(path: str | Path) -> Iterator[Frame]:
    source = Path(path)
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        yield from parse_bonds_lines(handle)


def parse_bonds_lines(lines: Iterable[str]) -> Iterator[Frame]:
    timestep: int | None = None
    expected_atom_count: int | None = None
    atoms: dict[int, AtomRecord] = {}
    output_cutoff: float | None = None
    source_line_start: int | None = None

    def finish_frame(line_number: int) -> Frame:
        if timestep is None:
            raise BondsParseError("Internal parser error: frame has no timestep")
        if expected_atom_count is not None and len(atoms) != expected_atom_count:
            raise BondsParseError(
                f"Timestep {timestep} declares {expected_atom_count} atoms but "
                f"contains {len(atoms)} data rows (near line {line_number})"
            )
        if not atoms and expected_atom_count != 0:
            raise BondsParseError(f"Timestep {timestep}: empty frame without an explicit zero atom count")
        for atom in atoms.values():
            # Check only this atom's neighbors, rather than traversing every
            # atom ID again for each row in the frame.
            missing = {neighbor for neighbor in atom.neighbors if neighbor not in atoms}
            if missing:
                raise BondsParseError(f"Timestep {timestep}: atom {atom.atom_id} references missing atoms {sorted(missing)}")
        return Frame(
            timestep=timestep,
            atoms=dict(sorted(atoms.items())),
            expected_atom_count=expected_atom_count,
            output_cutoff=output_cutoff,
            source_line_start=source_line_start,
            source_line_end=line_number,
        )

    last_line_number = 0
    for line_number, raw_line in enumerate(lines, start=1):
        last_line_number = line_number
        line = raw_line.strip()
        if not line:
            continue

        timestep_match = TIMESTEP_RE.match(line)
        if timestep_match:
            if timestep is not None:
                yield finish_frame(line_number - 1)
            timestep = int(timestep_match.group(1))
            expected_atom_count = None
            output_cutoff = None
            source_line_start = line_number
            atoms = {}
            continue

        atom_count_match = ATOM_COUNT_RE.match(line)
        if atom_count_match and timestep is not None:
            value = int(atom_count_match.group(1))
            if expected_atom_count is not None and expected_atom_count != value:
                raise BondsParseError(f"Conflicting atom counts at line {line_number}")
            expected_atom_count = value
            continue

        cutoff_match = OUTPUT_CUTOFF_RE.search(line) if line.startswith("#") else None
        if cutoff_match and timestep is not None:
            try:
                value = float(cutoff_match.group(1))
                if not math.isfinite(value) or value < 0:
                    raise ValueError("invalid output cutoff")
                if output_cutoff is not None and output_cutoff != value:
                    raise ValueError("conflicting output cutoffs")
                output_cutoff = value
            except ValueError as exc:
                raise BondsParseError(f"Invalid output cutoff at line {line_number}") from exc
            continue

        if line.startswith("#"):
            continue
        if timestep is None:
            raise BondsParseError(
                f"Atom data appears before the first '# Timestep' header at line {line_number}"
            )

        words = line.split()
        try:
            atom_id = int(words[0])
            atom_type = int(words[1])
            neighbor_count = int(words[2])
            if atom_id <= 0 or atom_type <= 0 or neighbor_count < 0:
                raise ValueError("invalid atom ID, type or neighbor count")
            minimum_fields = 7 + 2 * neighbor_count
            if len(words) < minimum_fields:
                raise ValueError(
                    f"expected at least {minimum_fields} fields, found {len(words)}"
                )
            neighbors = tuple(int(value) for value in words[3 : 3 + neighbor_count])
            molecule_id = int(words[3 + neighbor_count])
            bond_orders = tuple(
                float(value) for value in words[4 + neighbor_count : 4 + 2 * neighbor_count]
            )
            abo = float(words[4 + 2 * neighbor_count])
            lone_pairs = float(words[5 + 2 * neighbor_count])
            charge = float(words[6 + 2 * neighbor_count])
            if any(n <= 0 or n == atom_id for n in neighbors) or len(set(neighbors)) != len(neighbors):
                raise ValueError("invalid, self or duplicate neighbor")
            if not all(math.isfinite(v) for v in (*bond_orders, abo, lone_pairs, charge)):
                raise ValueError("non-finite atom value")
            if any(v < 0 for v in bond_orders):
                raise ValueError("negative bond order")
        except (IndexError, ValueError) as exc:
            raise BondsParseError(f"Malformed atom row at line {line_number}: {line}") from exc

        if atom_id in atoms:
            raise BondsParseError(
                f"Atom id {atom_id} occurs twice in timestep {timestep} at line {line_number}"
            )
        atoms[atom_id] = AtomRecord(
            atom_id=atom_id,
            atom_type=atom_type,
            neighbors=neighbors,
            molecule_id=molecule_id,
            bond_orders=bond_orders,
            abo=abo,
            lone_pairs=lone_pairs,
            charge=charge,
        )

    if timestep is not None:
        yield finish_frame(last_line_number)
