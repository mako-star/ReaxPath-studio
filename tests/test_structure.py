from __future__ import annotations

import csv
from contextlib import closing, redirect_stderr
import io
import json
import sqlite3
import sys
import tempfile
import unittest
from collections import Counter
from dataclasses import asdict, replace
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem, rdBase

from reaxnet.connectivity import ConnectivityConfig, run_connectivity
from reaxnet.identity import identify_structure_candidate
from reaxnet.structure import MolecularGraph, StructureConfig, infer_structure
from reaxnet.structures import main
from test_reaxnet import _frame, _trajectory


def graph_from_smiles(smiles):
    with rdBase.BlockLogs():
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return MolecularGraph(tuple(a.GetSymbol() for a in mol.GetAtoms()),
                          tuple((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()))


def graph_frame(graph, timestep=0, reverse=False, order=0.6):
    ids = list(range(1, len(graph.symbols) + 1))
    if reverse:
        ids.reverse()
    types = {ids[i]: {'C': 1, 'H': 2, 'O': 3}[symbol] for i, symbol in enumerate(graph.symbols)}
    edges = {(ids[a], ids[b]): order for a, b in graph.edges}
    return _frame(timestep, edges, atom_count=len(ids), atom_types=types)


class StructureSearchTests(unittest.TestCase):
    def test_known_CHO_structures_preserve_atoms_and_identify(self):
        for smiles in ('C', '[CH3]', 'C=O', 'CC', 'C=C', 'C#C', 'CCO', 'COC', '[H]', 'O', '[OH]', 'O=C=O', '[H][H]'):
            with self.subTest(smiles=smiles):
                graph = graph_from_smiles(smiles)
                result = infer_structure(graph, StructureConfig((0,), 1))
                self.assertEqual(result.identity_status, 'supported', asdict(result))
                self.assertTrue(result.search_complete)
                self.assertEqual({c.identifiers.canonical_smiles for c in result.candidates}, {smiles})
                for candidate in result.candidates:
                    self.assertEqual(candidate.identifiers.validation_status, 'passed')
                    with rdBase.BlockLogs():
                        restored = Chem.AddHs(Chem.MolFromSmiles(candidate.identifiers.canonical_smiles))
                    self.assertEqual(Counter(a.GetSymbol() for a in restored.GetAtoms()), Counter(graph.symbols))

    def test_hydrogen_edges_are_never_promoted(self):
        graph = graph_from_smiles('C=O')
        result = infer_structure(graph, StructureConfig((0,), 0))
        candidate = result.candidates[0]
        for edge, order in zip(graph.edges, candidate.bond_orders):
            self.assertEqual(order, 1 if 'H' in [graph.symbols[a] for a in edge] else 2)

    def test_radical_budget_changes_identity_not_hydrogen_count(self):
        methyl = graph_from_smiles('[CH3]')
        closed = infer_structure(methyl, StructureConfig((0,), 0))
        self.assertEqual(closed.identity_status, 'unresolved')
        self.assertEqual(closed.reason, 'no_candidate_under_declared_model')
        self.assertEqual(len(closed.candidates), 0)
        result = infer_structure(graph_from_smiles('C=O'), StructureConfig((0,), 2))
        self.assertEqual(result.identity_status, 'ambiguous')
        self.assertEqual({sum(c.radical_electrons) for c in result.candidates}, {0, 2})
        self.assertTrue(result.search_complete)

    def test_charge_unknown_methyl_retains_three_explanations(self):
        graph = graph_from_smiles('[CH3]')
        result = infer_structure(graph, StructureConfig((-1, 0, 1), 1))
        self.assertEqual(result.identity_status, 'ambiguous')
        self.assertEqual({c.identifiers.canonical_smiles for c in result.candidates}, {'[CH3]', '[CH3+]', '[CH3-]'})
        for charge in (-1, 1):
            charged = infer_structure(graph, StructureConfig((charge,), 0))
            self.assertEqual(charged.identity_status, 'supported')
            self.assertEqual(sum(charged.candidates[0].formal_charges), charge)

    def test_charge_separation_is_an_explicit_policy(self):
        graph = graph_from_smiles('C=O')
        neutral = infer_structure(graph, StructureConfig((0,), 0))
        separated = infer_structure(graph, StructureConfig((0,), 0, allow_charge_separation=True))
        self.assertEqual(neutral.identity_status, 'supported')
        self.assertEqual(separated.identity_status, 'ambiguous')
        self.assertIn('[CH2+][O-]', {c.identifiers.canonical_smiles for c in separated.candidates})
        self.assertNotEqual(neutral.assumption_id, separated.assumption_id)

    def test_kekule_candidates_merge_only_after_standardization(self):
        result = infer_structure(graph_from_smiles('c1ccccc1'), StructureConfig((0,), 0))
        self.assertEqual(result.identity_status, 'supported')
        self.assertEqual(len(result.candidates), 2)
        self.assertNotEqual(result.candidates[0].bond_orders, result.candidates[1].bond_orders)
        self.assertEqual({c.identifiers.canonical_smiles for c in result.candidates}, {'c1ccccc1'})

    def test_candidate_and_node_limits_cannot_claim_uniqueness(self):
        graph = graph_from_smiles('C=O')
        result = infer_structure(graph, StructureConfig((0,), 2, max_candidates=1))
        self.assertEqual(len(result.candidates), 1)
        self.assertEqual(result.identity_status, 'unresolved')
        self.assertFalse(result.search_complete)
        self.assertEqual(result.reason, 'max_candidates_exceeded')
        result = infer_structure(graph, StructureConfig((0,), 0, max_search_nodes=1))
        self.assertFalse(result.search_complete)
        self.assertEqual(result.search_nodes, 1)
        # Reaching the storage count exactly is not itself an incomplete search.
        result = infer_structure(graph_from_smiles('C'), StructureConfig((0,), 0, max_candidates=1))
        self.assertTrue(result.search_complete)
        self.assertEqual(result.identity_status, 'supported')

    def test_unsupported_large_and_impossible_graphs_remain_unresolved(self):
        for graph, config, reason in [
            (MolecularGraph(('N',), ()), StructureConfig((0,), 1), 'unsupported_elements'),
            (graph_from_smiles('CC'), StructureConfig((0,), 0, max_atoms=2), 'max_atoms_exceeded'),
            (graph_from_smiles('CCC'), StructureConfig((0,), 0, max_variable_edges=1), 'max_variable_edges_exceeded'),
            (MolecularGraph(('C',) + ('H',)*5, tuple((0,i) for i in range(1,6))), StructureConfig((0,), 1), 'no_candidate_under_declared_model'),
        ]:
            with self.subTest(reason=reason):
                result = infer_structure(graph, config)
                self.assertEqual(result.identity_status, 'unresolved')
                self.assertEqual(result.reason, reason)

    def test_missing_rdkit_preserves_candidates(self):
        with patch.dict(sys.modules, {'rdkit': None}):
            result = infer_structure(graph_from_smiles('C'), StructureConfig((0,), 0))
        self.assertEqual(result.identity_status, 'unresolved')
        self.assertEqual(len(result.candidates), 1)
        self.assertEqual(result.candidates[0].identifiers.smiles_status, 'unavailable')

    def test_inchi_failure_does_not_erase_supported_smiles(self):
        with patch('rdkit.Chem.inchi.MolToInchi', side_effect=RuntimeError('forced InChI failure')):
            result = infer_structure(graph_from_smiles('C'), StructureConfig((0,), 0))
        self.assertEqual(result.identity_status, 'supported')
        identifier = result.candidates[0].identifiers
        self.assertEqual(identifier.canonical_smiles, 'C')
        self.assertEqual(identifier.inchi_status, 'failed')
        self.assertIn('forced InChI failure', identifier.message)

    def test_invalid_candidate_is_not_repaired_with_added_hydrogens(self):
        graph = graph_from_smiles('[CH3]')
        identity = identify_structure_candidate(graph.symbols, graph.edges, (1,)*3, (0,)*4, (0,)*4)
        self.assertEqual(identity.validation_status, 'failed')
        self.assertEqual(identity.canonical_smiles, '')

    def test_identifier_disagreement_is_diagnostic(self):
        original = identify_structure_candidate
        calls = []
        def inconsistent(*args):
            value = original(*args)
            calls.append(value)
            return replace(value, inchikey='TEST-KEY-' + str(len(calls)))
        with patch('reaxnet.structure.identify_structure_candidate', side_effect=inconsistent):
            result = infer_structure(graph_from_smiles('c1ccccc1'), StructureConfig((0,), 0))
        self.assertEqual(result.identity_status, 'unresolved')
        self.assertEqual(result.reason, 'identifier_disagreement_for_same_smiles')
        self.assertEqual(len(result.candidates), 2)

    def test_invalid_configuration_is_rejected(self):
        for kwargs in ({'total_charges': ()}, {'total_charges': (0.5,)}, {'max_radicals': -1},
                       {'max_search_nodes': 0}, {'max_candidates': 0}, {'max_atoms': 500}, {'cache_size': -1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                StructureConfig(**{'total_charges': (0,), 'max_radicals': 1, **kwargs})


class StructurePipelineTests(unittest.TestCase):
    def run_fixture(self, root, text, structure_config):
        source, out = root / 'input.bonds', root / 'run'
        source.write_text(text)
        manifest = run_connectivity(source, {1:'C',2:'H',3:'O'}, out,
                                    ConnectivityConfig(), structure_config=structure_config)
        return out, manifest

    def test_renumbering_merges_species_and_preserves_local_atom_mapping(self):
        with tempfile.TemporaryDirectory() as d:
            graph = graph_from_smiles('CCO')
            # Two runs allow arbitrary ID/type relabeling without simulating transmutation.
            ids = []
            for index, reverse in enumerate((False, True)):
                root = Path(d)/str(index)
                root.mkdir()
                out, manifest = self.run_fixture(root, graph_frame(graph, reverse=reverse), StructureConfig((0,), 0))
                with (out/'species_registry.csv').open() as handle:
                    registry = list(csv.DictReader(handle))
                ids.append(registry[0]['species_id'])
                candidate_set = json.loads((out/'identity_candidates.jsonl').read_text())
                instance = json.loads((out/'instance_identities.jsonl').read_text())
                frame = json.loads((out/'frames.jsonl').read_text())
                types = {a['atom_id']: a['atom_type'] for a in frame['atoms']}
                self.assertEqual(candidate_set['local_graph']['symbols'], [{1:'C',2:'H',3:'O'}[types[a]] for a in instance['atom_ids']])
                self.assertEqual(manifest['structure']['summary']['instance_status_counts'], {'supported':1})
            self.assertEqual(ids[0], ids[1])

    def test_isomers_do_not_merge_and_low_BO_still_supports_double_bond(self):
        ids = []
        for smiles in ('CCO', 'COC', 'C=O'):
            with tempfile.TemporaryDirectory() as d:
                out, _ = self.run_fixture(Path(d), graph_frame(graph_from_smiles(smiles), order=0.6), StructureConfig((0,),0))
                with (out/'species_registry.csv').open() as handle:
                    row = next(csv.DictReader(handle))
                self.assertEqual(row['canonical_smiles'], smiles)
                ids.append(row['species_id'])
        self.assertNotEqual(ids[0], ids[1])

    def test_ambiguous_instance_is_counted_once_and_indexed(self):
        with tempfile.TemporaryDirectory() as d:
            out, manifest = self.run_fixture(Path(d), graph_frame(graph_from_smiles('[CH3]')), StructureConfig((-1,0,1),1))
            with (out/'species_timeseries.csv').open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows),1)
            self.assertEqual(rows[0]['identity_status'], 'ambiguous')
            self.assertEqual(rows[0]['count'], '1')
            evidence = json.loads((out/'identity_candidates.jsonl').read_text())
            self.assertEqual(len(evidence['candidates']),3)
            with closing(sqlite3.connect(out/'structure_index.sqlite')) as db:
                self.assertEqual(db.execute('SELECT count(*) FROM instances JOIN species USING(species_id) JOIN candidate_sets USING(candidate_set_id)').fetchone()[0],1)
            self.assertEqual(manifest['status'],'complete')

    def test_cache_is_bounded_and_repeated_instances_still_count(self):
        graph = graph_from_smiles('C')
        text = _trajectory(*(graph_frame(graph, timestep=i*10, order=0.6 + i*0.01) for i in range(20)))
        with tempfile.TemporaryDirectory() as d:
            out, manifest = self.run_fixture(Path(d),text,StructureConfig((0,),0,cache_size=1))
            summary=manifest['structure']['summary']
            self.assertEqual(summary['peak_cache_entries'],1)
            self.assertEqual(summary['inference_calls'],1)
            self.assertEqual(summary['cache_hits'],19)
            self.assertEqual(summary['instance_status_counts'],{'supported':20})
            self.assertEqual(len((out/'identity_candidates.jsonl').read_text().splitlines()),1)
            self.assertEqual(len((out/'instance_identities.jsonl').read_text().splitlines()),20)

    def test_lru_evicts_without_changing_species_or_counts(self):
        methane_edges = {(1, i): 0.9 for i in range(2,6)}
        methyl_edges = {(1, i): 0.9 for i in range(2,5)}
        text = _trajectory(*(_frame(i*10, methane_edges if i % 2 == 0 else methyl_edges) for i in range(6)))
        results = []
        for capacity in (0,1,10):
            with tempfile.TemporaryDirectory() as d:
                out, manifest = self.run_fixture(Path(d), text, StructureConfig((0,),1,cache_size=capacity))
                summary = manifest['structure']['summary']
                self.assertLessEqual(summary['peak_cache_entries'], capacity)
                with (out/'species_registry.csv').open() as h:
                    registry = {r['species_id']:r['canonical_smiles'] for r in csv.DictReader(h)}
                with (out/'species_timeseries.csv').open() as h:
                    results.append([(r['frame_id'],registry[r['species_id']],int(r['count'])) for r in csv.DictReader(h)])
        self.assertEqual(sorted(results[0]),sorted(results[1]))
        self.assertEqual(sorted(results[0]),sorted(results[2]))

    def test_failed_pipeline_marks_manifest_incomplete(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            graph=graph_from_smiles('C')
            with patch('reaxnet.species.infer_structure',side_effect=RuntimeError('forced failure')):
                with self.assertRaisesRegex(RuntimeError,'forced failure'):
                    self.run_fixture(root,graph_frame(graph),StructureConfig((0,),0))
            manifest=json.loads((root/'run'/'run_manifest.json').read_text())
            self.assertEqual(manifest['status'],'incomplete')
            self.assertIn('forced failure',manifest['error'])

    def test_cli_requires_declared_electronic_conditions(self):
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as ctx:
            main(['--bonds','none','--elements','1:C','--output-dir','none'])
        self.assertEqual(ctx.exception.code,2)


if __name__ == '__main__':
    unittest.main()
