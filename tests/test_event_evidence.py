from __future__ import annotations

import csv
from contextlib import closing
import json
import sqlite3
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from reaxnet.analysis import AnalyzerConfig, analyze_frames
from reaxnet.parser import parse_bonds_lines
from reaxnet.streaming import analyze_file_streaming
from test_reaxnet import _frame, _trajectory, CH4_BONDS, CH3_H_BONDS, ELEMENTS


class EventEvidenceTests(unittest.TestCase):
    def analyze(self, text, **kwargs):
        return analyze_frames(parse_bonds_lines(text.splitlines()), ELEMENTS,
                              AnalyzerConfig(generate_identifiers=False, **kwargs))

    def test_full_atoms_interval_and_conservation(self):
        result = self.analyze(_trajectory(_frame(100, CH4_BONDS), _frame(110, CH3_H_BONDS)), timestep_scale=0.001)
        evidence = result.events[0].evidence
        self.assertEqual(evidence['time_start_ps'], 0.1)
        self.assertEqual(evidence['time_end_ps'], 0.11)
        self.assertEqual(evidence['conservation'], {'atom_ids': True, 'elements': True})
        self.assertEqual(len(evidence['reactants'][0]['atoms']), 5)
        self.assertEqual(sorted(a['atom_id'] for m in evidence['products'] for a in m['atoms']), [1,2,3,4,5])
        self.assertEqual(evidence['reactants'][0]['connections'], [[1,2],[1,3],[1,4],[1,5]])
        self.assertNotEqual(evidence['reactants'][0]['instance_id'], evidence['products'][0]['instance_id'])

    def test_independent_simultaneous_events_are_separate(self):
        types = {1:1,2:2,3:1,4:2}
        text = _trajectory(_frame(0,{(1,2):0.4,(3,4):0.4},4,types),_frame(10,{},4,types))
        result = self.analyze(text)
        self.assertEqual(len(result.events), 2)
        self.assertEqual({e.involved_atom_ids for e in result.events}, {(1,2),(3,4)})
        self.assertTrue(all(e.evidence['conservation']['atom_ids'] for e in result.events))

    def test_type_change_and_time_gap_are_rejected(self):
        cases = [
            _trajectory(_frame(0,CH4_BONDS),_frame(10,CH4_BONDS,atom_types={2:1})),
            _trajectory(_frame(0,CH4_BONDS),_frame(10,CH4_BONDS),_frame(30,CH3_H_BONDS)),
        ]
        for text in cases:
            with self.subTest(text=text), self.assertRaises(ValueError): self.analyze(text)
        with self.assertRaises(ValueError):
            self.analyze(_trajectory(_frame(0,CH4_BONDS),_frame(20,CH4_BONDS)), expected_interval=10)

    def test_isomerization_same_formula_keeps_structure_difference(self):
        types={1:1,2:1,3:1,4:1}
        text=_trajectory(_frame(0,{(1,2):1,(2,3):1,(3,4):1},4,types),
                         _frame(10,{(1,2):1,(1,3):1,(1,4):1},4,types))
        event=self.analyze(text).events[0]
        self.assertEqual(event.reactant_formulas,event.product_formulas)
        self.assertNotEqual(event.reactant_keys,event.product_keys)
        self.assertEqual(event.event_type,'isomerization')

    def test_identical_structure_with_different_ids_merges(self):
        types={1:1,2:2,3:2,4:1,5:2,6:2}
        result=self.analyze(_frame(0,{(1,2):0.4,(1,3):0.4,(4,5):0.4,(4,6):0.4},6,types))
        self.assertEqual(len(result.registry),1)

    def test_renumbering_with_changed_element_layout_uses_exact_isomorphism(self):
        types={1:1,2:1,3:2,4:2,5:2,6:1,7:2,8:1}
        edges={(1,2):1,(1,3):1,(2,4):1,(6,8):1,(6,5):1,(8,7):1}
        result=self.analyze(_frame(0,edges,8,types))
        self.assertEqual(len(result.registry),1)
        first,second=result.states[0].molecules
        self.assertNotEqual(first.atom_elements,second.atom_elements)
        self.assertEqual(first.structure_id,'C2H2-001')
        self.assertEqual(first.structure_id,second.structure_id)

    def test_fingerprint_collision_does_not_merge_different_graphs(self):
        types={i:1 for i in range(1,9)}
        edges={(1,2):1,(2,3):1,(3,4):1,(5,6):1,(5,7):1,(5,8):1}
        with patch('reaxnet.chemistry._graph_fingerprint', return_value='collision'):
            result=self.analyze(_frame(0,edges,8,types))
        self.assertEqual(len(result.registry),2)

    def test_index_payload_matches_csv_final_recrossing_status(self):
        text=_trajectory(_frame(0,CH4_BONDS),_frame(10,CH3_H_BONDS),_frame(20,CH4_BONDS))
        with tempfile.TemporaryDirectory() as d:
            source=Path(d)/'input.bonds';source.write_text(text)
            result=analyze_file_streaming(source,ELEMENTS,Path(d)/'out',
                       AnalyzerConfig(generate_identifiers=False,timestep_scale=0.001,recrossing_window=0.03))
            with closing(sqlite3.connect(result.paths['event_index'])) as db:
                rows=db.execute('select event_id,accepted,filter_reason,payload from events order by event_id').fetchall()
                self.assertEqual(len(rows),2)
                self.assertGreater(db.execute('select count(*) from event_species').fetchone()[0],0)
            with result.paths['reaction_events'].open() as f: csv_rows=list(csv.DictReader(f))
            for row,csv_row in zip(rows,csv_rows):
                payload=json.loads(row[3])
                self.assertEqual(row[0],csv_row['event_id'])
                self.assertFalse(payload['accepted'])
                self.assertEqual(payload['filter_reason'],csv_row['filter_reason'])
                self.assertEqual(payload['reaction_identifiers'],csv_row['reaction'])
                self.assertIn('-001', payload['reaction'])
                self.assertTrue(payload['conservation']['elements'])

    def test_structure_labels_distinguish_isomers_and_count_stoichiometry(self):
        from reaxnet.outputs import event_to_payload
        types={i:1 for i in range(1,5)}
        result=self.analyze(_trajectory(_frame(0,{(1,2):1,(2,3):1,(3,4):1},4,types),
                        _frame(10,{(1,2):1,(1,3):1,(1,4):1},4,types)))
        payload=event_to_payload(result.events[0],result.registry)
        self.assertEqual(payload['reaction'],'C4-001 -> C4-002')
        exchange_types={1:1,2:2,3:2,4:1,5:2,6:2}
        result=self.analyze(_trajectory(
            _frame(0,{(1,2):1,(1,3):1,(4,5):1,(4,6):1},6,exchange_types),
            _frame(10,{(1,5):1,(1,3):1,(4,2):1,(4,6):1},6,exchange_types)))
        self.assertEqual(event_to_payload(result.events[0],result.registry)['reaction'],
                         '2 CH2-001 -> 2 CH2-001')

    def test_identity_failure_does_not_rewrite_topology(self):
        from reaxnet.identity import identify_molecule
        from rdkit import Chem
        from reaxnet.models import Bond
        bonds={(1,2):Bond(1,2,3.0,3.0)}
        with patch.object(Chem,'SanitizeMol',side_effect=ValueError('invalid bond interpretation')) as sanitize:
            result=identify_molecule((1,2),{1:'C',2:'H'},bonds)
        self.assertEqual(result.status,'failed')
        self.assertEqual(result.smiles,'')
        self.assertEqual(sanitize.call_count,1)
        self.assertEqual(set(bonds),{(1,2)})

    def test_identifier_cache_separates_bond_interpretations(self):
        from reaxnet.identity import IdentityResult
        text=_trajectory(_frame(0,{(1,2):0.4},2),_frame(10,{(1,2):2.0},2))
        with patch('reaxnet.chemistry.identify_molecule',side_effect=[
            IdentityResult(smiles='first'),IdentityResult(smiles='second')]) as identify:
            result=analyze_frames(parse_bonds_lines(text.splitlines()),ELEMENTS,AnalyzerConfig())
        self.assertEqual(identify.call_count,2)
        self.assertEqual(result.states[0].molecules[0].smiles,'first')
        self.assertEqual(result.states[1].molecules[0].smiles,'second')
        self.assertEqual(result.states[0].molecules[0].species_key,result.states[1].molecules[0].species_key)
        self.assertEqual(len(result.events),0)

    def test_exchange_closure_keeps_entire_unchanged_fragments(self):
        types={1:1,2:2,3:2,4:1,5:2,6:2}
        first={(1,2):1,(1,3):1,(4,5):1,(4,6):1}
        second={(1,5):1,(1,3):1,(4,2):1,(4,6):1}
        event=self.analyze(_trajectory(_frame(0,first,6,types),_frame(10,second,6,types))).events[0]
        self.assertEqual(len(event.evidence['reactants']),2)
        self.assertEqual(len(event.evidence['products']),2)
        self.assertEqual(event.involved_atom_ids,(1,2,3,4,5,6))
        self.assertEqual(event.evidence['conservation'],{'atom_ids':True,'elements':True})

    def test_exported_structure_ids_and_observed_interval(self):
        with tempfile.TemporaryDirectory() as d:
            source=Path(d)/'input.bonds'
            source.write_text(_trajectory(_frame(0,CH4_BONDS),_frame(10,CH3_H_BONDS)))
            result=analyze_file_streaming(source,ELEMENTS,Path(d)/'out',AnalyzerConfig(generate_identifiers=False))
            self.assertEqual(result.summary['observed_interval'],10)
            self.assertEqual(result.network['metadata']['observed_interval'],10)
            with result.paths['molecule_registry'].open() as f:
                registry={r['species_key']:r['structure_id'] for r in csv.DictReader(f)}
            self.assertEqual(set(registry.values()),{'CH4-001','CH3-001','H-001'})
            with result.paths['species_timeseries'].open() as f:
                rows=list(csv.DictReader(f))
            self.assertTrue(rows)
            self.assertTrue(all(r['structure_id']==registry[r['species_key']] for r in rows))

    def test_streaming_rejects_missing_frame(self):
        with tempfile.TemporaryDirectory() as d:
            source=Path(d)/'input.bonds';source.write_text(_trajectory(_frame(0,CH4_BONDS),_frame(10,CH4_BONDS),_frame(30,CH3_H_BONDS)))
            with self.assertRaises(ValueError):
                analyze_file_streaming(source,ELEMENTS,Path(d)/'out',AnalyzerConfig(generate_identifiers=False))

if __name__ == '__main__': unittest.main()
