"""Streaming element-pair counts and explicit local-topology group candidates.

File adjacency is authoritative.  Functional-group counts are motif matches,
not inferred bond orders, formal charges, or confirmed chemical identities.
"""
from __future__ import annotations

import csv
from collections import Counter
from collections.abc import Callable, Mapping
from contextlib import ExitStack
from dataclasses import dataclass
from itertools import combinations_with_replacement, islice
import json
import math
from pathlib import Path
from time import monotonic
from typing import Any

from reaxnet.connectivity import ConnectivityConfig, ConnectivityFrame, ConnectivityTracker
from reaxnet.models import Frame
from reaxnet.parser import canonical_symbol, parse_bonds_file


SCHEMA_VERSION = "trajectory-counts-1"
RULE_VERSION = "cho-local-topology-candidates-1"
BOND_CATEGORIES = ("C-C", "C-H", "C-O", "O-H", "H-H", "O-O")
GROUP_KEYS = ("hydroxyl", "ether", "aldehyde", "ketone", "carboxyl", "ester")
GROUP_LABELS = {
    "hydroxyl": "羟基候选", "ether": "醚候选", "aldehyde": "醛基候选",
    "ketone": "酮型羰基候选", "carboxyl": "羧基候选", "ester": "酯基候选",
}
COUNT_FILENAMES = {
    "bond_counts": "bond_counts.csv", "bond_changes": "bond_changes.csv",
    "functional_groups": "functional_groups.csv",
    "functional_group_molecules": "functional_group_molecules.csv",
    "count_metadata": "count_analysis.json",
}

RULES = {
    "basis": "validated_file_adjacency_only",
    "chemical_identity": "unconfirmed_local_topology_candidates",
    "bond_orders": "not_used; no secondary BO cutoff or single/double/triple classification",
    "count_unit": "one motif per distinct center atom and group; carriers are distinct connected components",
    "hydroxyl": "O has exactly C,H neighbors; H has degree 1; exclude the OH assigned to a carboxyl candidate. Alcohol/phenol not separated.",
    "ether": "O has exactly two C neighbors; exclude the bridging O assigned to an ester candidate.",
    "aldehyde": "C has exactly three neighbors: one terminal O and two C/H atoms, with at least one H. H atoms have degree 1.",
    "ketone": "C has exactly three neighbors: one terminal O and two C atoms. This is a ketone-like candidate, not all carbonyl groups.",
    "carboxyl": "C has exactly three neighbors: two O and one C/H; one O is terminal and the other has exactly C(center),H neighbors.",
    "ester": "C has exactly three neighbors: two O and one C/H; one O is terminal and the other has exactly C(center),C(other) neighbors.",
    "terminal_oxygen": "O has exactly one neighbor, the candidate carbon center; no C=O bond order is asserted.",
    "priority": "carboxyl and ester before hydroxyl and ether; an acid OH is not counted as ordinary hydroxyl and ester O is not ordinary ether",
    "local_exclusion": "skip motifs involving unsupported atoms or atoms with observed degree above H:1,C:4,O:2; keep all such edges in bond counts",
    "small_molecules": "water, isolated OH, CO and CO2 do not match these organic-group rules",
    "zero": "zero means no match under these rules, not proof that a chemical group is absent",
    "limits": "Topology alone cannot resolve bond multiplicity, resonance, charge or radical state. Carbonyl-like motifs can have alternative electronic interpretations. Do not add category counts as a chemically exhaustive total.",
}


def bond_pair_label(first: str, second: str) -> str:
    pair = frozenset((first, second))
    if pair == {"H", "O"}:
        return "O-H"
    return "-".join(sorted((first, second)))


def bond_categories(elements: Mapping[int, str]) -> tuple[str, ...]:
    symbols = sorted(set(elements.values()))
    extra = {bond_pair_label(a, b) for a, b in combinations_with_replacement(symbols, 2)}
    return BOND_CATEGORIES + tuple(sorted(extra - set(BOND_CATEGORIES)))


@dataclass(frozen=True)
class FunctionalGroupSite:
    group: str
    center_atom_id: int
    atom_ids: tuple[int, ...]
    molecule_index: int


def functional_group_sites(
    frame: Frame, observed: ConnectivityFrame, elements: Mapping[int, str],
) -> tuple[tuple[FunctionalGroupSite, ...], dict[str, int]]:
    """Match strict local motifs without changing or relabeling any connection."""
    symbols = {aid: elements[atom.atom_type] for aid, atom in frame.atoms.items()}
    adjacency = {aid: set() for aid in frame.atoms}
    membership = {}
    for index, molecule in enumerate(observed.instances):
        membership.update((aid, index) for aid in molecule.atom_ids)
        for first, second in molecule.connections:
            adjacency[first].add(second)
            adjacency[second].add(first)
    degree_limits = {"H": 1, "C": 4, "O": 2}
    unsupported = {aid for aid, symbol in symbols.items() if symbol not in degree_limits}
    overcoordinated = {aid for aid, neighbors in adjacency.items()
                       if symbols[aid] in degree_limits and len(neighbors) > degree_limits[symbols[aid]]}
    excluded = unsupported | overcoordinated
    sites = []
    claimed_oxygen = set()
    seen = set()

    def add(group, center, atoms):
        atoms = set(atoms)
        if atoms & excluded or (group, center) in seen:
            return
        seen.add((group, center))
        sites.append(FunctionalGroupSite(group, center, tuple(sorted(atoms)), membership[center]))

    def neighbors_of_type(atom_id, symbol):
        return [nid for nid in adjacency[atom_id] if symbols[nid] == symbol]

    # Larger patterns claim their OH/bridging O before generic O-centered motifs.
    for center, symbol in symbols.items():
        if symbol != "C" or center in excluded or len(adjacency[center]) != 3:
            continue
        neighbors = adjacency[center]
        if neighbors & excluded:
            continue
        oxygens = neighbors_of_type(center, "O")
        if len(oxygens) != 2:
            continue
        remainder = neighbors - set(oxygens)
        if len(remainder) != 1 or symbols[next(iter(remainder))] not in {"C", "H"}:
            continue
        terminals = [oid for oid in oxygens if adjacency[oid] == {center}]
        if len(terminals) != 1:
            continue
        terminal = terminals[0]
        bridge = next(oid for oid in oxygens if oid != terminal)
        others = adjacency[bridge] - {center}
        if len(adjacency[bridge]) != 2 or len(others) != 1 or others & excluded:
            continue
        external = next(iter(others))
        group = {"H": "carboxyl", "C": "ester"}.get(symbols[external])
        if group:
            add(group, center, neighbors | {center, external})
            claimed_oxygen.add(bridge)

    for center, symbol in symbols.items():
        neighbors = adjacency[center]
        if center in excluded or neighbors & excluded:
            continue
        if symbol == "O" and len(neighbors) == 2 and center not in claimed_oxygen:
            types = Counter(symbols[nid] for nid in neighbors)
            if types == {"C": 1, "H": 1}:
                add("hydroxyl", center, neighbors | {center})
            elif types == {"C": 2}:
                add("ether", center, neighbors | {center})
        elif symbol == "C" and len(neighbors) == 3:
            oxygens = neighbors_of_type(center, "O")
            if len(oxygens) != 1 or adjacency[oxygens[0]] != {center}:
                continue
            other_types = Counter(symbols[nid] for nid in neighbors if nid != oxygens[0])
            if set(other_types) <= {"C", "H"} and other_types["H"] >= 1:
                add("aldehyde", center, neighbors | {center})
            elif other_types == {"C": 2}:
                add("ketone", center, neighbors | {center})

    matched_oxygen = {aid for site in sites for aid in site.atom_ids if symbols[aid] == "O"}
    oxygen_atoms = {aid for aid, symbol in symbols.items() if symbol == "O"}
    coverage = {
        "atom_observations": len(symbols),
        "unsupported_atom_observations": len(unsupported),
        "overcoordinated_atom_observations": len(overcoordinated),
        "matched_oxygen_atom_observations": len(matched_oxygen),
        "unmatched_oxygen_atom_observations": len(oxygen_atoms - matched_oxygen),
        "frames_with_unsupported_elements": int(bool(unsupported)),
        "frames_with_overcoordination": int(bool(overcoordinated)),
    }
    return tuple(sites), coverage


def _source_stat(path: Path) -> dict[str, int]:
    stat = path.stat()
    return {"size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


class CountWriter:
    """Write only the active frame, reusing caller-validated connectivity."""

    def __init__(self, output_dir, element_mapping, timestep_scale, *, source_path=None,
                 expected_interval=None):
        if isinstance(timestep_scale, bool) or not math.isfinite(timestep_scale) or timestep_scale <= 0:
            raise ValueError("timestep_scale must be finite and positive")
        ConnectivityConfig(expected_interval=expected_interval)
        if not element_mapping:
            raise ValueError("an explicit element mapping is required")
        self.elements = {key: canonical_symbol(value) for key, value in element_mapping.items()}
        self.scale = timestep_scale
        self.expected_interval = expected_interval
        self.interval = expected_interval
        self.source = None if source_path is None else Path(source_path).resolve()
        self.initial_stat = None if self.source is None else _source_stat(self.source)
        self.output_dir = Path(output_dir)
        self.paths = {key: self.output_dir / name for key, name in COUNT_FILENAMES.items()}
        self.categories = bond_categories(self.elements)
        self.previous_edges = None
        self.previous_timestep = None
        self.first_timestep = None
        self.frame_count = 0
        self.coverage = Counter()
        self._stack = None
        self._writers = {}
        self.metadata: dict[str, Any] = {}

    @property
    def summary(self):
        return {"frame_count": self.frame_count, "first_timestep": self.first_timestep,
                "last_timestep": self.previous_timestep, "observed_interval": self.interval}

    def _save_metadata(self, status, error=None):
        self.metadata = {
            "schema_version": SCHEMA_VERSION, "rule_version": RULE_VERSION, "status": status,
            **self.summary,
            "time_start_ps": None if self.first_timestep is None else self.first_timestep * self.scale,
            "time_end_ps": None if self.previous_timestep is None else self.previous_timestep * self.scale,
            "timestep_scale_ps": self.scale, "expected_interval": self.expected_interval,
            "elements": self.elements,
            "bond_categories": [{"key": key, "label": key.replace("-", "–")} for key in self.categories],
            "functional_group_categories": [{"key": key, "label": GROUP_LABELS[key]} for key in GROUP_KEYS],
            "rules": RULES, "coverage": {
                **dict(self.coverage), "supported_group_elements": ["C", "H", "O"],
                "unmatched_oxygen_meaning": "not matched by these organic-group rules; includes water, gases and unclassified local environments",
            },
            "source": None if self.source is None else str(self.source),
            "source_stat": self.initial_stat, "error": error,
        }
        target = self.paths["count_metadata"]
        temporary = target.with_suffix(".tmp")
        temporary.write_text(json.dumps(self.metadata, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
        temporary.replace(target)

    def __enter__(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)
        if any(path.exists() for path in self.paths.values()):
            raise FileExistsError("Count outputs already exist; use a new output directory")
        self._stack = ExitStack()
        base = ["timestep", "time_ps"]
        fields = {
            "bond_counts": base + ["total"] + list(self.categories),
            "bond_changes": base + ["previous_timestep", "previous_time_ps", "total_formed", "total_broken"] +
                [name for key in self.categories for name in (f"{key}_formed", f"{key}_broken")],
            "functional_groups": base + list(GROUP_KEYS),
            "functional_group_molecules": base + list(GROUP_KEYS),
        }
        try:
            for key, columns in fields.items():
                handle = self._stack.enter_context(self.paths[key].open("x", newline="", encoding="utf-8"))
                writer = csv.DictWriter(handle, fieldnames=columns)
                writer.writeheader()
                self._writers[key] = writer
            self._save_metadata("writing")
        except BaseException:
            self._stack.close()
            raise
        return self

    def write(self, frame: Frame, observed: ConnectivityFrame):
        if self._stack is None or not self._writers:
            raise RuntimeError("CountWriter must be opened before writing")
        if frame.timestep != observed.timestep:
            raise ValueError("Connectivity frame does not match timestep")
        if self.previous_timestep is not None:
            interval = frame.timestep - self.previous_timestep
            if interval <= 0 or (self.interval is not None and interval != self.interval):
                raise ValueError("Count frames must have increasing, consistent timesteps")
            if self.interval is None:
                self.interval = interval
        symbols = {aid: self.elements[atom.atom_type] for aid, atom in frame.atoms.items()}
        edges = {edge for molecule in observed.instances for edge in molecule.connections}

        def histogram(selected):
            counts = Counter(bond_pair_label(symbols[a], symbols[b]) for a, b in selected)
            return {key: counts[key] for key in self.categories}

        base = {"timestep": frame.timestep, "time_ps": frame.timestep * self.scale}
        stock = {**base, "total": len(edges), **histogram(edges)}
        delta = {**base, "previous_timestep": self.previous_timestep,
                 "previous_time_ps": None if self.previous_timestep is None else self.previous_timestep * self.scale}
        formed = None if self.previous_edges is None else edges - self.previous_edges
        broken = None if self.previous_edges is None else self.previous_edges - edges
        delta.update(total_formed=None if formed is None else len(formed),
                     total_broken=None if broken is None else len(broken))
        for suffix, selection in (("formed", formed), ("broken", broken)):
            values = None if selection is None else histogram(selection)
            delta.update({f"{key}_{suffix}": None if values is None else values[key] for key in self.categories})
        sites, coverage = functional_group_sites(frame, observed, self.elements)
        group_counts = Counter(site.group for site in sites)
        carriers = {key: set() for key in GROUP_KEYS}
        for site in sites:
            carriers[site.group].add(site.molecule_index)
        self._writers["bond_counts"].writerow(stock)
        self._writers["bond_changes"].writerow(delta)
        self._writers["functional_groups"].writerow({**base, **{key: group_counts[key] for key in GROUP_KEYS}})
        self._writers["functional_group_molecules"].writerow({**base, **{key: len(carriers[key]) for key in GROUP_KEYS}})
        self.coverage.update(coverage)
        self.first_timestep = frame.timestep if self.first_timestep is None else self.first_timestep
        self.previous_timestep = frame.timestep
        self.previous_edges = edges
        self.frame_count += 1

    def __exit__(self, exc_type, exc, traceback):
        self._stack.close()
        self._writers.clear()
        if exc_type is not None:
            self._save_metadata("failed", f"{exc_type.__name__}: {exc}")
            return False
        try:
            if self.source is not None and _source_stat(self.source) != self.initial_stat:
                raise ValueError("Source file changed during count analysis")
            self._save_metadata("completed" if self.frame_count else "incomplete",
                                None if self.frame_count else "No frames written")
        except BaseException as error:
            self._save_metadata("failed", f"{type(error).__name__}: {error}")
            raise
        return False


@dataclass(frozen=True)
class CountAnalysisResult:
    paths: dict[str, Path]
    summary: dict[str, Any]
    metadata: dict[str, Any]


def analyze_counts_file(source, element_mapping, output_dir, *, timestep_scale=1.0,
                        expected_interval=None, max_frames=None,
                        progress_callback: Callable[[dict], None] | None = None):
    """Read a bounded source once, without species isomorphism or event analysis."""
    if max_frames is not None and (isinstance(max_frames, bool) or not isinstance(max_frames, int) or max_frames <= 0):
        raise ValueError("max_frames must be a positive integer")
    source = Path(source).resolve()
    writer = CountWriter(output_dir, element_mapping, timestep_scale, source_path=source,
                         expected_interval=expected_interval)
    tracker = ConnectivityTracker(element_mapping, ConnectivityConfig(expected_interval=expected_interval))
    parser = parse_bonds_file(source)
    last_progress = 0.0
    with writer:
        try:
            frames = parser if max_frames is None else islice(parser, max_frames)
            for frame in frames:
                observed = tracker.process(frame)
                writer.write(frame, observed)
                now = monotonic()
                if progress_callback is not None and now - last_progress >= 0.5:
                    progress_callback({"frames_processed": writer.frame_count,
                                       "timestep": frame.timestep, "max_frames": max_frames})
                    last_progress = now
            if not writer.frame_count:
                raise ValueError("The bonds trajectory contains no timestep frames")
            if progress_callback is not None:
                progress_callback({"frames_processed": writer.frame_count,
                                   "timestep": writer.previous_timestep, "max_frames": max_frames})
        finally:
            parser.close()
    return CountAnalysisResult(writer.paths, writer.summary, writer.metadata)
