"""Disk-backed species registration and per-instance structure evidence."""
from __future__ import annotations

import csv
import json
import sqlite3
from collections import Counter, OrderedDict
from dataclasses import asdict
from pathlib import Path

from reaxnet.structure import MolecularGraph, StructureConfig, digest, infer_structure

REGISTRY_FIELDS = ('species_id', 'formula', 'canonical_smiles', 'inchi', 'inchikey',
                   'formal_charge', 'radical_electrons', 'identity_status', 'smiles_status',
                   'inchi_status', 'assumption_id', 'representative_candidate_set_id', 'message')
TIMESERIES_FIELDS = ('frame_id', 'timestep', 'time_ps', 'species_id', 'formula', 'count',
                     'identity_status', 'view')


class StructureWriter:
    """Bounded inference LRU; growing registry/instance indices are kept on disk."""
    def __init__(self, target: Path, config: StructureConfig):
        self.target, self.config = Path(target), config
        self.assumption = config.assumption_id
        self.cache = OrderedDict()
        self.status_counts = Counter()
        self.cache_hits = 0
        self.inference_calls = 0
        self.peak_cache_entries = 0
        self.handles = []
        self.db = None

    def __enter__(self):
        try:
            self.db = sqlite3.connect(self.target / 'structure_index.sqlite')
            self.db.executescript('''
                CREATE TABLE species (species_id TEXT PRIMARY KEY, payload TEXT NOT NULL);
                CREATE TABLE candidate_sets (candidate_set_id TEXT PRIMARY KEY, evidence_line INTEGER NOT NULL);
                CREATE TABLE instances (instance_id TEXT PRIMARY KEY, frame_id TEXT NOT NULL,
                    species_id TEXT NOT NULL, candidate_set_id TEXT NOT NULL);
                CREATE INDEX instances_by_frame ON instances(frame_id);
                CREATE INDEX instances_by_species ON instances(species_id);
            ''')
            def open_output(name):
                handle = (self.target / name).open('w', encoding='utf-8', newline='')
                self.handles.append(handle)
                return handle
            self.candidates_file = open_output('identity_candidates.jsonl')
            self.instances_file = open_output('instance_identities.jsonl')
            self.timeseries = csv.DictWriter(open_output('species_timeseries.csv'), fieldnames=TIMESERIES_FIELDS)
            self.timeseries.writeheader()
            self.candidate_lines = 0
            return self
        except BaseException:
            self._close()
            raise

    def _close(self):
        for handle in self.handles:
            handle.close()
        if self.db is not None:
            self.db.close()

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            if exc_type is None:
                self.db.commit()
                with (self.target / 'species_registry.csv').open('w', encoding='utf-8', newline='') as handle:
                    writer = csv.DictWriter(handle, fieldnames=REGISTRY_FIELDS)
                    writer.writeheader()
                    for (payload,) in self.db.execute('SELECT payload FROM species ORDER BY species_id'):
                        writer.writerow(json.loads(payload))
        finally:
            self._close()
        return False

    def _infer(self, graph):
        # BO, partial q and geometry are deliberately not inputs to model v1.
        # The exact ordered graph prevents unsafe reuse between different local
        # atom layouts. Isomorphic layouts can miss the cache, but never mis-map.
        key = (graph, self.assumption)
        if key in self.cache:
            self.cache_hits += 1
            result = self.cache.pop(key)
            self.cache[key] = result
            return result
        result = infer_structure(graph, self.config)
        self.inference_calls += 1
        if self.config.cache_size:
            self.cache[key] = result
            while len(self.cache) > self.config.cache_size:
                self.cache.popitem(last=False)
            self.peak_cache_entries = max(self.peak_cache_entries, len(self.cache))
        return result

    def process(self, raw_frame, connectivity_frame, elements):
        symbols = {a: elements[record.atom_type] for a, record in raw_frame.atoms.items()}
        counts = Counter()
        rows = {}
        for instance in connectivity_frame.instances:
            graph = MolecularGraph.from_instance(instance, symbols)
            result = self._infer(graph)
            candidate_set_id = 'CS' + digest({'graph': asdict(graph), 'assumption': self.assumption})
            if self.db.execute('SELECT 1 FROM candidate_sets WHERE candidate_set_id=?', (candidate_set_id,)).fetchone() is None:
                self.candidate_lines += 1
                payload = {'candidate_set_id': candidate_set_id, 'local_graph': asdict(graph),
                           'index_convention': 'zero_based_positions_in_instance_atom_ids', **asdict(result)}
                self.candidates_file.write(json.dumps(payload, separators=(',', ':')) + '\n')
                self.db.execute('INSERT INTO candidate_sets VALUES (?,?)', (candidate_set_id, self.candidate_lines))
            if result.identity_status == 'supported':
                candidate = next((c for c in result.candidates if c.identifiers.inchi_status == 'ok'), result.candidates[0])
                identifier = candidate.identifiers
                identity_key = ['supported', self.assumption, identifier.canonical_smiles]
                fields = {'canonical_smiles': identifier.canonical_smiles, 'inchi': identifier.inchi,
                          'inchikey': identifier.inchikey, 'formal_charge': sum(candidate.formal_charges),
                          'radical_electrons': sum(candidate.radical_electrons),
                          'smiles_status': identifier.smiles_status, 'inchi_status': identifier.inchi_status}
            else:
                # Conservative graph buckets; never merge an uncertain graph
                # with a supported structure or expand its candidates as counts.
                identity_key = [result.identity_status, self.assumption, asdict(graph), result.reason]
                fields = {'canonical_smiles': '', 'inchi': '', 'inchikey': '', 'formal_charge': '',
                          'radical_electrons': '', 'smiles_status': 'candidate_specific', 'inchi_status': 'candidate_specific'}
            species_id = 'S' + digest(identity_key)
            row = {'species_id': species_id, 'formula': instance.formula, **fields,
                   'identity_status': result.identity_status, 'assumption_id': self.assumption,
                   'representative_candidate_set_id': candidate_set_id, 'message': result.reason}
            self.db.execute('INSERT OR IGNORE INTO species VALUES (?,?)', (species_id, json.dumps(row)))
            self.db.execute('INSERT INTO instances VALUES (?,?,?,?)',
                            (instance.instance_id, instance.frame_id, species_id, candidate_set_id))
            self.instances_file.write(json.dumps({
                'instance_id': instance.instance_id, 'frame_id': instance.frame_id, 'species_id': species_id,
                'candidate_set_id': candidate_set_id, 'atom_ids': instance.atom_ids,
                'identity_status': result.identity_status, 'search_complete': result.search_complete,
                'assumption_id': self.assumption, 'reason': result.reason,
            }, separators=(',', ':')) + '\n')
            counts[species_id] += 1
            rows[species_id] = row
            self.status_counts[result.identity_status] += 1
        for species_id, count in sorted(counts.items()):
            row = rows[species_id]
            self.timeseries.writerow({'frame_id': connectivity_frame.frame_id, 'timestep': raw_frame.timestep,
                                      'time_ps': connectivity_frame.time_ps, 'species_id': species_id,
                                      'formula': row['formula'], 'count': count,
                                      'identity_status': row['identity_status'], 'view': 'raw_instances'})
        self.db.commit()

    def summary(self):
        return {'instance_status_counts': dict(self.status_counts),
                'species_count': self.db.execute('SELECT count(*) FROM species').fetchone()[0],
                'candidate_set_count': self.candidate_lines, 'cache_hits': self.cache_hits,
                'inference_calls': self.inference_calls, 'peak_cache_entries': self.peak_cache_entries}
