from __future__ import annotations

import json
import io
from contextlib import redirect_stderr
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from reaxnet.connectivity import ConnectivityConfig, ConnectivityTracker, build_connectivity_parser, run_connectivity
from reaxnet.parser import BondsParseError, parse_bonds_lines
from test_reaxnet import _frame, _trajectory


def frame(timestep=0, order=None):
    edges = {} if order is None else {(1, 2): order}
    return next(parse_bonds_lines(_frame(timestep, edges, atom_count=2).splitlines()))


def tracker(**kwargs):
    return ConnectivityTracker({1: 'C', 2: 'H'}, ConnectivityConfig(**kwargs))


class InputValidationTests(unittest.TestCase):
    def test_preserves_cutoff_and_source_lines(self):
        f = frame(order=0.9)
        self.assertEqual(f.output_cutoff, 0.3)
        self.assertEqual(f.source_line_start, 1)
        self.assertGreater(f.source_line_end, 8)

    def test_malformed_rows_are_rejected(self):
        header = '# Timestep 0\n# Number of particles 1\n'
        rows = [
            '1 1 -1 0 0 0 0',
            '1 1 0 0 nan 0 0',
            '1 1 1 1 0 1 1 0 0',
            '1 1 2 2 2 0 1 1 2 0 0',
            '1 1 1 2 0 1 1 0 0',
            '0 1 0 0 0 0 0',
            '1 1 1 2 0 -1 -1 0 0',
        ]
        for row in rows:
            with self.subTest(row=row), self.assertRaises(BondsParseError):
                list(parse_bonds_lines((header + row).splitlines()))

    def test_truncated_and_duplicate_atoms_are_rejected(self):
        for text in ['# Timestep 0', '# Timestep 0\n# Number of particles 2\n1 1 0 0 0 0 0',
                     '# Timestep 0\n1 1 0 0 0 0 0\n1 1 0 0 0 0 0']:
            with self.subTest(text=text), self.assertRaises(BondsParseError):
                list(parse_bonds_lines(text.splitlines()))


class ConnectivityTests(unittest.TestCase):
    def test_every_listed_edge_is_kept_regardless_of_bo(self):
        for order in (0.0, 0.299, 0.3, 0.35, 0.499, 0.5, 0.501, 0.9):
            with self.subTest(order=order):
                result = tracker().process(frame(order=order))
                self.assertTrue(result.observations[0].connected)
                self.assertEqual(result.observations[0].reason, 'listed_in_frame')
                self.assertEqual(result.instances[0].connections, ((1, 2),))
                self.assertEqual(len(result.instances), 1)
                self.assertIsNone(result.time_ps)

    def test_connectivity_tracks_presence_without_hysteresis(self):
        t = tracker()
        for i, order in enumerate((0.4, 0.8, 0.35, None, 0.3)):
            with self.subTest(order=order):
                result = t.process(frame(i * 10, order))
                self.assertEqual(t.connected, set() if order is None else {(1, 2)})
                self.assertEqual(len(result.instances), 2 if order is None else 1)

    def test_missing_record_retains_censored_evidence(self):
        t = tracker()
        t.process(frame(0, 0.8))
        result = t.process(frame(10))
        observation = result.observations[0]
        self.assertIsNone(observation.merged_order)
        self.assertEqual(observation.directed_orders, ())
        self.assertFalse(observation.connected)
        self.assertEqual(len(result.instances), 2)
        explicit_zero = tracker().process(frame(order=0.0)).observations[0]
        self.assertEqual(explicit_zero.merged_order, 0.0)

    def test_output_cutoff_is_optional_metadata(self):
        t = tracker()
        result = t.process(replace(frame(order=0.3), output_cutoff=None))
        self.assertIsNone(result.output_cutoff)
        self.assertEqual(len(result.instances), 1)
        absent = t.process(replace(frame(10), output_cutoff=None)).observations[0]
        self.assertFalse(absent.connected)
        self.assertIsNone(absent.merged_order)
        self.assertEqual(absent.reason, 'not_listed_in_frame')
        self.assertEqual(tracker(output_cutoff=0.8).process(
            replace(frame(order=0.3), output_cutoff=None)).instances[0].connections, ((1, 2),))
        with self.assertRaisesRegex(ValueError, 'disagrees'):
            tracker(output_cutoff=0.4).process(frame())

    def test_conflicting_and_one_sided_observations(self):
        original = frame(order=0.6)
        conflict = replace(original, atoms={**original.atoms, 2: replace(original.atoms[2], bond_orders=(0.7,))})
        with self.assertRaisesRegex(ValueError, 'conflicting BO'):
            tracker().process(conflict)
        result = tracker(symmetry_tolerance=0.1).process(conflict)
        self.assertAlmostEqual(result.observations[0].merged_order, 0.65)
        self.assertEqual(result.observations[0].directed_orders, ((1, 0.6), (2, 0.7)))
        one_sided = replace(original, atoms={**original.atoms, 2: replace(original.atoms[2], neighbors=(), bond_orders=())})
        with self.assertRaisesRegex(ValueError, 'one directed'):
            tracker().process(one_sided)

    def test_frame_local_identity(self):
        t = tracker(timestep_scale=0.0002)
        a = t.process(frame(0, 0.6))
        b = t.process(frame(1000, 0.6))
        self.assertEqual([m.formula for m in a.instances], ['CH'])
        self.assertNotEqual(a.instances[0].instance_id, b.instances[0].instance_id)
        self.assertEqual(a.instances[0].atom_ids, b.instances[0].atom_ids)
        self.assertEqual(b.time_ps, 0.2)

    def test_connectivity_does_not_delete_fifth_carbon_neighbor(self):
        raw = next(parse_bonds_lines(_frame(0, {(1, a): 0.9 for a in range(2, 7)}, atom_count=6).splitlines()))
        result = tracker().process(raw)
        self.assertEqual(len(result.instances), 1)
        self.assertEqual(result.instances[0].formula, 'CH5')
        self.assertEqual(len(result.instances[0].connections), 5)
        # This is a graph formula, not a claim that CH5 is a supported species.

    def test_gap_type_changes_and_missing_mapping_do_not_advance_state(self):
        t = tracker(expected_interval=10)
        first = t.process(frame(0, 0.6))
        invalid = [frame(0, 0.6), frame(-10, 0.6), frame(20, 0.6),
                   replace(frame(10), atoms={1: frame().atoms[1]}, expected_atom_count=1),
                   replace(frame(10), atoms={**frame().atoms, 2: replace(frame().atoms[2], atom_type=1)})]
        for f in invalid:
            with self.subTest(timestep=f.timestep), self.assertRaises(ValueError):
                t.process(f)
            self.assertEqual(t.frame_count, 1)
        self.assertEqual(t.process(frame(10, 0.6)).frame_id, 'F00000001')
        with self.assertRaisesRegex(ValueError, 'Missing element'):
            ConnectivityTracker({1: 'C'}, ConnectivityConfig()).process(frame())

    def test_invalid_numeric_config(self):
        for number in [float('nan'), float('inf'), -1]:
            with self.subTest(number=number), self.assertRaises(ValueError):
                tracker(symmetry_tolerance=number)
        with self.assertRaises(ValueError):
            tracker(timestep_scale=float('nan'))


class RunContractTests(unittest.TestCase):
    def test_complete_run_is_auditable_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            source, out = root / 'input.bonds', root / 'run'
            source.write_text(_trajectory(_frame(0, {}, atom_count=2), _frame(10, {(1, 2): 0.9}, atom_count=2)))
            config = ConnectivityConfig()
            result = run_connectivity(source, {1: 'C', 2: 'H'}, out, config)
            self.assertEqual(result['status'], 'complete')
            self.assertEqual(result['algorithm_version'], 'file-topology-2')
            self.assertEqual(result['policy']['connectivity'], 'file_neighbor_list')
            self.assertNotIn('thresholds', result['config'])
            self.assertEqual(result['frame_count'], 2)
            instances = [json.loads(s) for s in (out / 'molecule_instances.jsonl').read_text().splitlines()]
            frames = [json.loads(s) for s in (out / 'frames.jsonl').read_text().splitlines()]
            self.assertEqual(len(instances), 3)
            self.assertEqual(frames[1]['first_instance_line'], 3)
            self.assertEqual(instances[2]['connections'], [[1, 2]])
            self.assertEqual(frames[1]['atoms'][0]['bond_orders'], [0.9])
            self.assertEqual(instances[2]['frame_id'], frames[1]['frame_id'])
            with self.assertRaises(FileExistsError):
                run_connectivity(source, {1: 'C', 2: 'H'}, out, config)

    def test_removed_threshold_options_are_rejected(self):
        from reaxnet.cli import build_parser
        required = ['--bonds', 'input.bonds', '--elements', '1:C', '--output-dir', 'run']
        for parser, options in (
            (build_connectivity_parser(), ['--formation', '0.5']),
            (build_connectivity_parser(), ['--breaking', '0.35']),
            (build_connectivity_parser(), ['--pair', 'C-H:0.5']),
            (build_parser(), ['--bond-cutoff', '0.5']),
            (build_parser(), ['--valence-overrides', 'C:4']),
        ):
            with self.subTest(options=options), redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                parser.parse_args(required + options)
            self.assertEqual(error.exception.code, 2)

    def test_partial_and_failed_runs_never_claim_complete(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            source = root / 'input.bonds'
            source.write_text(_trajectory(_frame(0, {}, atom_count=2), _frame(20, {}, atom_count=2)))
            config = ConnectivityConfig(expected_interval=10)
            self.assertEqual(run_connectivity(source, {1: 'C', 2: 'H'}, root / 'partial', config, max_frames=1)['status'], 'partial')
            with self.assertRaisesRegex(ValueError, 'missing frame'):
                run_connectivity(source, {1: 'C', 2: 'H'}, root / 'failed', config)
            manifest = json.loads((root / 'failed' / 'run_manifest.json').read_text())
            self.assertEqual(manifest['status'], 'incomplete')
            self.assertEqual(manifest['frame_count'], 1)
            self.assertIn('error', manifest)


if __name__ == '__main__':
    unittest.main()
