"""PRD stage 1: observed connectivity and frame-local molecule instances.

No chemical bond types, formal charges or species identities are inferred here.
"""
from __future__ import annotations

import argparse
from contextlib import ExitStack
import json
import math
import platform
from time import perf_counter
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping
from uuid import uuid4

from reaxnet.models import BondKey, Frame
from reaxnet.parser import canonical_symbol, parse_bonds_file, parse_element_mapping


@dataclass(frozen=True)
class ConnectivityConfig:
    symmetry_tolerance: float = 0.0
    output_cutoff: float | None = None
    expected_interval: int | None = None
    timestep_scale: float | None = None

    def __post_init__(self):
        if not math.isfinite(self.symmetry_tolerance) or self.symmetry_tolerance < 0:
            raise ValueError("symmetry_tolerance must be finite and nonnegative")
        if self.output_cutoff is not None and (
            not math.isfinite(self.output_cutoff) or self.output_cutoff < 0
        ):
            raise ValueError("output_cutoff must be finite and nonnegative")
        if self.timestep_scale is not None and (
            not math.isfinite(self.timestep_scale) or self.timestep_scale <= 0
        ):
            raise ValueError("timestep_scale must be finite and positive")
        if self.expected_interval is not None and (
            isinstance(self.expected_interval, bool)
            or not isinstance(self.expected_interval, int) or self.expected_interval <= 0
        ):
            raise ValueError("expected_interval must be a positive integer")


@dataclass(frozen=True)
class PairObservation:
    atom_ids: BondKey
    # Directed values are retained; None means not listed, never BO=0.
    directed_orders: tuple[tuple[int, float], ...]
    merged_order: float | None
    connected: bool
    reason: str


@dataclass(frozen=True)
class MoleculeInstance:
    instance_id: str
    frame_id: str
    atom_ids: tuple[int, ...]
    formula: str
    connections: tuple[BondKey, ...]


@dataclass(frozen=True)
class ConnectivityFrame:
    frame_id: str
    timestep: int
    time_ps: float | None
    source_line_start: int | None
    source_line_end: int | None
    output_cutoff: float | None
    observations: tuple[PairObservation, ...]
    instances: tuple[MoleculeInstance, ...]


def _components(atom_ids, edges):
    adjacency = {a: [] for a in atom_ids}
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    seen = set()
    for start in sorted(adjacency):
        if start in seen:
            continue
        seen.add(start)
        stack, component = [start], []
        while stack:
            a = stack.pop()
            component.append(a)
            for b in adjacency[a]:
                if b not in seen:
                    seen.add(b)
                    stack.append(b)
        yield tuple(sorted(component))


def _formula(atoms, symbols):
    counts = Counter(symbols[a] for a in atoms)
    order = ([s for s in ('C', 'H') if s in counts] + sorted(counts.keys() - {'C', 'H'})
             if 'C' in counts else sorted(counts))
    return ''.join(s + (str(counts[s]) if counts[s] != 1 else '') for s in order)


class ConnectivityTracker:
    """Fail closed on invalid frames; commit history only after successful validation."""
    def __init__(self, element_mapping: Mapping[int, str], config: ConnectivityConfig):
        if not element_mapping:
            raise ValueError("an explicit element mapping is required")
        self.elements = {k: canonical_symbol(v) for k, v in element_mapping.items()}
        self.config = config
        self.previous_types = None
        self.previous_timestep = None
        self.interval = config.expected_interval
        self.connected: set[BondKey] = set()
        self.frame_count = 0

    def process(self, frame: Frame) -> ConnectivityFrame:
        if frame.expected_atom_count is None or len(frame.atoms) != frame.expected_atom_count:
            raise ValueError(f"Timestep {frame.timestep}: a matching declared atom count is required")
        types = {a: record.atom_type for a, record in frame.atoms.items()}
        if self.previous_types is not None and types != self.previous_types:
            raise ValueError(f"Timestep {frame.timestep}: atom IDs or types changed; start a declared segment")
        interval = self.interval
        if self.previous_timestep is not None:
            delta = frame.timestep - self.previous_timestep
            if delta <= 0 or (interval is not None and delta != interval):
                raise ValueError(f"Timestep {frame.timestep}: duplicate, reversed or missing frame; expected interval {interval}")
            interval = delta if interval is None else interval
        missing = set(types.values()) - self.elements.keys()
        if missing:
            raise ValueError(f"Missing element mapping for types {sorted(missing)}")
        symbols = {a: self.elements[t] for a, t in types.items()}
        cutoff = frame.output_cutoff
        if self.config.output_cutoff is not None:
            if cutoff is not None and cutoff != self.config.output_cutoff:
                raise ValueError("configured output cutoff disagrees with frame header")
            cutoff = self.config.output_cutoff
        if cutoff is not None and (not math.isfinite(cutoff) or cutoff < 0):
            raise ValueError("output cutoff must be finite and nonnegative when known")
        # The neighbor list is authoritative. The writer cutoff is metadata only.
        directed = {}
        for atom_id, atom in frame.atoms.items():
            if atom.atom_id != atom_id or atom_id <= 0 or atom.atom_type <= 0:
                raise ValueError("invalid atom ID or type")
            if len(atom.neighbors) != len(atom.bond_orders) or len(set(atom.neighbors)) != len(atom.neighbors):
                raise ValueError("invalid or duplicate neighbor records")
            for neighbor, order in zip(atom.neighbors, atom.bond_orders):
                if neighbor == atom_id or neighbor not in frame.atoms or not math.isfinite(order) or order < 0:
                    raise ValueError("invalid neighbor or bond order")
                key = tuple(sorted((atom_id, neighbor)))
                directed.setdefault(key, []).append((atom_id, order))
        observations, connected = [], set()
        for key in sorted(directed.keys() | self.connected):
            values = tuple(sorted(directed.get(key, [])))
            if not values:
                order, accepted, reason = None, False, 'not_listed_in_frame'
            else:
                if len(values) != 2:
                    raise ValueError(f"Timestep {frame.timestep}: pair {key} has only one directed observation")
                difference = abs(values[0][1] - values[1][1])
                if difference > self.config.symmetry_tolerance + 1e-12:
                    raise ValueError(f"Timestep {frame.timestep}: pair {key} has conflicting BO values {values}")
                order = values[0][1] / 2 + values[1][1] / 2
                accepted, reason = True, 'listed_in_frame'
                if difference:
                    reason += ':mean_within_declared_tolerance'
            if accepted:
                connected.add(key)
            observations.append(PairObservation(key, values, order, accepted, reason))
        frame_id = f'F{self.frame_count:08d}'
        components = list(_components(frame.atoms, connected))
        membership = {a: i for i, component in enumerate(components) for a in component}
        component_edges = [[] for _ in components]
        for edge in sorted(connected):
            component_edges[membership[edge[0]]].append(edge)
        instances = tuple(MoleculeInstance(
            f'{frame_id}:M{i:08d}', frame_id, atoms, _formula(atoms, symbols), tuple(component_edges[i])
        ) for i, atoms in enumerate(components))
        result = ConnectivityFrame(
            frame_id, frame.timestep,
            None if self.config.timestep_scale is None else frame.timestep * self.config.timestep_scale,
            frame.source_line_start, frame.source_line_end, cutoff, tuple(observations), instances,
        )
        self.connected, self.previous_types = connected, types
        self.previous_timestep, self.interval = frame.timestep, interval
        self.frame_count += 1
        return result


def run_connectivity(source, elements, output_dir, config, *, max_frames=None, structure_config=None):
    """Write a new run directory. Partial and failed runs never claim completion."""
    started = perf_counter()
    if max_frames is not None and (isinstance(max_frames, bool) or not isinstance(max_frames, int) or max_frames <= 0):
        raise ValueError("max_frames must be a positive integer")
    source = Path(source).resolve()
    stat = source.stat()
    target = Path(output_dir).resolve()
    tracker = ConnectivityTracker(elements, config)
    target.mkdir(parents=True, exist_ok=False)
    manifest = {
        'schema_version': 'connectivity-2', 'algorithm_version': 'file-topology-2',
        'run_id': str(uuid4()), 'status': 'incomplete',
        'runtime': {'python': platform.python_version(), 'platform': platform.platform(), 'dependencies': 'stdlib_only'},
        'input': {'path': str(source), 'size_bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns,
                  'fingerprint_policy': 'path_size_mtime_only_no_content_hash'},
        'elements': tracker.elements,
        'config': asdict(config),
        'policy': {'connectivity': 'file_neighbor_list', 'output_cutoff': 'metadata_only', 'bidirectional_merge': 'mean_within_tolerance_else_error',
                   'one_sided_record': 'error', 'invalid_frame': 'stop_without_bridging',
                   'frame_interval': 'explicit_or_inferred_from_first_two_frames',
                   'chemical_identity': 'not_inferred', 'physical_time': 'ps' if config.timestep_scale else None},
        'max_frames': max_frames, 'frame_count': 0, 'instance_count': 0,
    }
    if structure_config is not None:
        from reaxnet.structure import MODEL_VERSION, NORMALIZATION_VERSION, rdkit_version
        manifest['structure'] = {'config': asdict(structure_config),
                                 'assumption_id': structure_config.assumption_id,
                                 'model_version': MODEL_VERSION, 'normalization_version': NORMALIZATION_VERSION,
                                 'rdkit_version': rdkit_version(), 'zero_count_policy': 'omit_zero_rows',
                                 'stereochemistry': 'unspecified_no_coordinates',
                                 'bo_and_partial_charge_usage': 'evidence_only_not_used_for_inference',
                                 'charge_policy': 'per_fragment_explicit_not_system_charge',
                                 'unknown_identity_policy': 'separate_conservative_ordered_graph_buckets'}
        manifest['schema_version'] = 'structures-2'
        manifest['policy']['chemical_identity'] = 'bounded_CHO_candidate_search'
        manifest['runtime']['dependencies'] = 'stdlib_and_optional_rdkit'
    manifest_path = target / 'run_manifest.json'
    def save_manifest():
        temporary = target / 'run_manifest.tmp'
        temporary.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        temporary.replace(manifest_path)
    save_manifest()
    try:
        with ExitStack() as stack:
            frames_out = stack.enter_context((target / 'frames.jsonl').open('w', encoding='utf-8'))
            instances_out = stack.enter_context((target / 'molecule_instances.jsonl').open('w', encoding='utf-8'))
            structure_writer = None
            if structure_config is not None:
                from reaxnet.species import StructureWriter
                structure_writer = stack.enter_context(StructureWriter(target, structure_config))
            iterator = parse_bonds_file(source)
            try:
                for raw in iterator:
                    frame = tracker.process(raw)
                    frame_data = asdict(frame)
                    frame_data.pop('instances')
                    frame_data['atoms'] = [asdict(raw.atoms[a]) for a in sorted(raw.atoms)]
                    frame_data['instance_count'] = len(frame.instances)
                    frame_data['first_instance_line'] = manifest['instance_count'] + 1
                    frames_out.write(json.dumps(frame_data, separators=(',', ':'), allow_nan=False) + '\n')
                    for instance in frame.instances:
                        instances_out.write(json.dumps(asdict(instance), separators=(',', ':')) + '\n')
                    if structure_writer is not None:
                        structure_writer.process(raw, frame, tracker.elements)
                    manifest['frame_count'] += 1
                    manifest['instance_count'] += len(frame.instances)
                    if max_frames is not None and manifest['frame_count'] >= max_frames:
                        manifest['status'] = 'partial'
                        break
                else:
                    if not manifest['frame_count']:
                        raise ValueError('no frames found')
                    manifest['status'] = 'complete'
            finally:
                iterator.close()
                if structure_writer is not None:
                    manifest['structure']['summary'] = structure_writer.summary()
        after = source.stat()
        if (stat.st_size, stat.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError('input changed during analysis')
        manifest['observed_interval'] = tracker.interval
        manifest['elapsed_seconds'] = perf_counter() - started
        manifest['data_output_bytes'] = sum(p.stat().st_size for p in target.iterdir() if p.is_file() and p.name != 'run_manifest.json')
        save_manifest()
    except BaseException as exc:
        manifest['status'] = 'incomplete'
        manifest['elapsed_seconds'] = perf_counter() - started
        manifest['error'] = f'{type(exc).__name__}: {exc}'
        save_manifest()
        raise
    return manifest


def build_connectivity_parser():
    parser = argparse.ArgumentParser(description='PRD stage 1: connectivity and molecule instances')
    parser.add_argument('--bonds', type=Path, required=True)
    parser.add_argument('--elements', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--symmetry-tolerance', type=float, default=0.0)
    parser.add_argument('--output-cutoff', type=float, help='Writer cutoff metadata only; does not filter neighbors')
    parser.add_argument('--expected-interval', type=int)
    parser.add_argument('--timestep-scale', type=float, help='Explicit ps per simulation step; omitted means no physical time')
    parser.add_argument('--max-frames', type=int, help='Sample only; manifest is marked partial')
    return parser


def connectivity_config_from_args(args):
    return ConnectivityConfig(symmetry_tolerance=args.symmetry_tolerance,
                              output_cutoff=args.output_cutoff,
                              expected_interval=args.expected_interval,
                              timestep_scale=args.timestep_scale)


def main(argv=None):
    parser = build_connectivity_parser()
    args = parser.parse_args(argv)
    try:
        config = connectivity_config_from_args(args)
        result = run_connectivity(args.bonds, parse_element_mapping(args.elements), args.output_dir,
                                  config, max_frames=args.max_frames)
    except (ValueError, OSError) as exc:
        parser.exit(2, f'connectivity: error: {exc}\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
