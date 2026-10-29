from __future__ import annotations

from collections import Counter
import csv
from dataclasses import replace
import json
from pathlib import Path
import tempfile
import unittest

from reaxnet.connectivity import ConnectivityConfig, ConnectivityTracker
from reaxnet.counts import CountWriter, analyze_counts_file, functional_group_sites
from reaxnet.models import AtomRecord, Frame
from test_reaxnet import _frame, _trajectory


ELEMENTS = {1: "C", 2: "H", 3: "O", 4: "N"}


def frame_from_graph(symbols, edges, timestep=0):
    reverse = {symbol: key for key, symbol in ELEMENTS.items()}
    adjacency = {aid: [] for aid in symbols}
    if not isinstance(edges, dict):
        edges = {pair: 1.0 for pair in edges}
    for (a, b), bo in edges.items():
        adjacency[a].append((b, bo))
        adjacency[b].append((a, bo))
    atoms = {aid: AtomRecord(aid, reverse[symbol], tuple(n for n, _ in adjacency[aid]), 0,
                            tuple(bo for _, bo in adjacency[aid]),
                            sum(bo for _, bo in adjacency[aid]), 0, 0)
             for aid, symbol in symbols.items()}
    return Frame(timestep, atoms, expected_atom_count=len(atoms))


def group_matches(frame):
    observed = ConnectivityTracker(ELEMENTS, ConnectivityConfig()).process(frame)
    return functional_group_sites(frame, observed, ELEMENTS)


def read_csv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class BondCountTests(unittest.TestCase):
    def test_undirected_counts_include_low_and_zero_bo_and_changes_conserve_stock(self):
        symbols = {1: "C", 2: "C", 3: "H", 4: "O"}
        first = frame_from_graph(symbols, {(1, 2): 0, (1, 3): 0.1, (2, 4): 2.7}, 10)
        second = frame_from_graph(symbols, {(1, 3): 0.1, (2, 4): 2.7, (3, 4): 0.01}, 20)
        tracker = ConnectivityTracker(ELEMENTS, ConnectivityConfig(expected_interval=10))
        with tempfile.TemporaryDirectory() as folder:
            with CountWriter(folder, ELEMENTS, 0.2) as writer:
                for frame in (first, second):
                    writer.write(frame, tracker.process(frame))
            counts = read_csv(writer.paths["bond_counts"])
            changes = read_csv(writer.paths["bond_changes"])
            self.assertEqual((counts[0]["C-C"], counts[0]["C-H"], counts[0]["C-O"]), ("1", "1", "1"))
            self.assertEqual(int(counts[0]["total"]), sum(len(a.neighbors) for a in first.atoms.values()) // 2)
            self.assertEqual((counts[0]["time_ps"], counts[1]["time_ps"]), ("2.0", "4.0"))
            self.assertEqual((changes[0]["total_formed"], changes[0]["C-C_broken"]), ("", ""))
            self.assertEqual((changes[1]["C-C_broken"], changes[1]["O-H_formed"]), ("1", "1"))
            self.assertEqual((changes[1]["total_formed"], changes[1]["total_broken"]), ("1", "1"))
            self.assertEqual(int(counts[1]["total"]), int(counts[0]["total"]) +
                             int(changes[1]["total_formed"]) - int(changes[1]["total_broken"]))
            self.assertEqual(writer.metadata["status"], "completed")
            self.assertEqual(writer.summary["observed_interval"], 10)

    def test_other_elements_get_counted_and_group_coverage_is_explicit(self):
        frame = frame_from_graph({1: "N", 2: "O", 3: "H"}, [(1, 2), (2, 3)])
        observed = ConnectivityTracker(ELEMENTS, ConnectivityConfig()).process(frame)
        with tempfile.TemporaryDirectory() as folder:
            with CountWriter(folder, ELEMENTS, 1) as writer:
                writer.write(frame, observed)
            row = read_csv(writer.paths["bond_counts"])[0]
            self.assertEqual((row["N-O"], row["O-H"], row["total"]), ("1", "1", "2"))
            self.assertEqual(writer.metadata["coverage"]["unsupported_atom_observations"], 1)
            self.assertEqual(writer.metadata["coverage"]["frames_with_unsupported_elements"], 1)
            self.assertIn("not proof", writer.metadata["rules"]["zero"])

    def test_empty_saved_frame_is_zero_stock_not_unknown_previous_change(self):
        frame = Frame(0, {}, expected_atom_count=0)
        observed = ConnectivityTracker(ELEMENTS, ConnectivityConfig()).process(frame)
        with tempfile.TemporaryDirectory() as folder:
            with CountWriter(folder, ELEMENTS, 1) as writer:
                writer.write(frame, observed)
            self.assertEqual(read_csv(writer.paths["bond_counts"])[0]["total"], "0")
            self.assertEqual(read_csv(writer.paths["functional_groups"])[0]["hydroxyl"], "0")
            self.assertEqual(read_csv(writer.paths["bond_changes"])[0]["total_formed"], "")
            self.assertEqual(writer.metadata["status"], "completed")
            self.assertIsNone(writer.metadata["source"])


class FunctionalGroupCountTests(unittest.TestCase):
    def test_two_hydroxyls_in_one_carrier_exclude_water_and_oh_radical(self):
        symbols = {1: "C", 2: "C", 3: "O", 4: "H", 5: "O", 6: "H",
                   7: "H", 8: "H", 9: "H", 10: "H", 11: "O", 12: "H", 13: "H",
                   14: "O", 15: "H"}
        edges = [(1, 2), (1, 3), (3, 4), (2, 5), (5, 6), (1, 7), (1, 8),
                 (2, 9), (2, 10), (11, 12), (11, 13), (14, 15)]
        frame = frame_from_graph(symbols, edges)
        sites, coverage = group_matches(frame)
        self.assertEqual(Counter(s.group for s in sites), {"hydroxyl": 2})
        self.assertEqual({s.center_atom_id for s in sites}, {3, 5})
        self.assertEqual(len({s.molecule_index for s in sites}), 1)
        self.assertEqual(coverage["unmatched_oxygen_atom_observations"], 2)
        with tempfile.TemporaryDirectory() as folder:
            with CountWriter(folder, ELEMENTS, 1) as writer:
                writer.write(frame, ConnectivityTracker(ELEMENTS, ConnectivityConfig()).process(frame))
            self.assertEqual(read_csv(writer.paths["functional_groups"])[0]["hydroxyl"], "2")
            self.assertEqual(read_csv(writer.paths["functional_group_molecules"])[0]["hydroxyl"], "1")

    def test_acid_and_ester_claim_their_oxygen_without_generic_double_counting(self):
        # Acetic acid plus methyl acetate, as separate complete components.
        acid_symbols = {1: "C", 2: "C", 3: "O", 4: "O", 5: "H", 6: "H", 7: "H", 8: "H"}
        acid_edges = [(1, 2), (2, 3), (2, 4), (4, 5), (1, 6), (1, 7), (1, 8)]
        ester_symbols = {20: "C", 21: "C", 22: "O", 23: "O", 24: "C",
                         25: "H", 26: "H", 27: "H", 28: "H", 29: "H", 30: "H"}
        ester_edges = [(20, 21), (21, 22), (21, 23), (23, 24), (20, 25), (20, 26),
                       (20, 27), (24, 28), (24, 29), (24, 30)]
        sites, _ = group_matches(frame_from_graph({**acid_symbols, **ester_symbols}, acid_edges + ester_edges))
        self.assertEqual(Counter(s.group for s in sites), {"carboxyl": 1, "ester": 1})
        self.assertEqual({s.center_atom_id for s in sites}, {2, 21})

    def test_formaldehyde_and_acetone_match_but_co_and_co2_do_not(self):
        symbols = {1: "C", 2: "O", 3: "H", 4: "H", 10: "C", 11: "O",
                   12: "C", 13: "C", 14: "H", 15: "H", 16: "H", 17: "H",
                   18: "H", 19: "H", 30: "C", 31: "O", 32: "O", 40: "C", 41: "O"}
        edges = [(1, 2), (1, 3), (1, 4), (10, 11), (10, 12), (10, 13),
                 (12, 14), (12, 15), (12, 16), (13, 17), (13, 18), (13, 19),
                 (30, 31), (30, 32), (40, 41)]
        sites, coverage = group_matches(frame_from_graph(symbols, edges))
        self.assertEqual(Counter(s.group for s in sites), {"aldehyde": 1, "ketone": 1})
        self.assertEqual(coverage["unmatched_oxygen_atom_observations"], 3)

    def test_ether_center_is_unique_even_with_symmetric_carbon_neighbors(self):
        symbols = {1: "O", 2: "C", 3: "C", **{i: "H" for i in range(4, 10)}}
        edges = [(1, 2), (1, 3), (2, 4), (2, 5), (2, 6), (3, 7), (3, 8), (3, 9)]
        sites, _ = group_matches(frame_from_graph(symbols, edges))
        self.assertEqual([(s.group, s.center_atom_id) for s in sites], [("ether", 1)])

    def test_overcoordinated_neighborhood_keeps_bonds_but_does_not_assert_group(self):
        symbols = {1: "C", 2: "O", 3: "H", **{i: "H" for i in range(4, 8)}}
        # C has five neighbors. Its O-H would look like an alcohol if only O
        # were checked; the local candidate is deliberately left unclassified.
        frame = frame_from_graph(symbols, [(1, 2), (2, 3), (1, 4), (1, 5), (1, 6), (1, 7)])
        sites, coverage = group_matches(frame)
        self.assertEqual(sites, ())
        self.assertEqual(coverage["overcoordinated_atom_observations"], 1)
        with tempfile.TemporaryDirectory() as folder:
            with CountWriter(folder, ELEMENTS, 1) as writer:
                writer.write(frame, ConnectivityTracker(ELEMENTS, ConnectivityConfig()).process(frame))
            self.assertEqual(read_csv(writer.paths["bond_counts"])[0]["total"], "6")
            self.assertEqual(read_csv(writer.paths["functional_groups"])[0]["hydroxyl"], "0")
            self.assertEqual(writer.metadata["coverage"]["frames_with_overcoordination"], 1)

    def test_four_coordinate_carbon_with_terminal_o_is_not_a_carbonyl_candidate(self):
        symbols = {1: "C", 2: "O", 3: "C", 4: "C", 5: "H"}
        sites, _ = group_matches(frame_from_graph(symbols, [(1, 2), (1, 3), (1, 4), (1, 5)]))
        self.assertEqual(sites, ())


class CountFileTests(unittest.TestCase):
    def test_bounded_stream_is_exact_and_does_not_need_reaction_analysis(self):
        mapping = {1: "C", 2: "H"}
        text = _trajectory(_frame(0, {(1, 2): 0.0}, 2), _frame(10, {}, 2))
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder) / "sample.bonds"
            source.write_text(text)
            progress = []
            result = analyze_counts_file(source, mapping, Path(folder)/"out", max_frames=1,
                                         timestep_scale=0.2, expected_interval=10,
                                         progress_callback=progress.append)
            self.assertEqual(result.summary["frame_count"], 1)
            self.assertEqual((result.metadata["first_timestep"], result.metadata["last_timestep"]), (0, 0))
            self.assertEqual(result.metadata["source_stat"]["size_bytes"], source.stat().st_size)
            self.assertEqual(len(read_csv(result.paths["bond_counts"])), 1)
            self.assertEqual(progress[-1]["frames_processed"], 1)
            with self.assertRaises(FileExistsError):
                analyze_counts_file(source, mapping, Path(folder)/"out")

    def test_invalid_limits_and_time_scale_fail_before_analysis(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder)/"sample.bonds"
            source.write_text(_frame(0, {}, 1))
            for limit in (0, -1, True, 1.5):
                with self.subTest(limit=limit), self.assertRaises(ValueError):
                    analyze_counts_file(source, ELEMENTS, Path(folder)/"out", max_frames=limit)
            for scale in (0, -1, True, float("nan"), float("inf")):
                with self.subTest(scale=scale), self.assertRaises(ValueError):
                    CountWriter(Path(folder)/"out", ELEMENTS, scale)

    def test_bad_interval_and_empty_input_never_mark_completed(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder)/"sample.bonds"
            source.write_text(_trajectory(_frame(0, {}, 1), _frame(20, {}, 1)))
            with self.assertRaises(ValueError):
                analyze_counts_file(source, ELEMENTS, Path(folder)/"gap", expected_interval=10)
            metadata = json.loads((Path(folder)/"gap/count_analysis.json").read_text())
            self.assertEqual((metadata["status"], metadata["frame_count"]), ("failed", 1))
            source.write_text("")
            with self.assertRaisesRegex(ValueError, "no timestep"):
                analyze_counts_file(source, ELEMENTS, Path(folder)/"empty")
            self.assertEqual(json.loads((Path(folder)/"empty/count_analysis.json").read_text())["status"], "failed")

    def test_asymmetric_adjacency_is_rejected_before_counting(self):
        frame = frame_from_graph({1: "C", 2: "H"}, [(1, 2)])
        asymmetric = replace(frame, atoms={1: frame.atoms[1],
            2: replace(frame.atoms[2], neighbors=(), bond_orders=())})
        with tempfile.TemporaryDirectory() as folder:
            with self.assertRaisesRegex(ValueError, "one directed"):
                with CountWriter(folder, ELEMENTS, 1) as writer:
                    writer.write(asymmetric, ConnectivityTracker(ELEMENTS, ConnectivityConfig()).process(asymmetric))
            self.assertEqual(writer.metadata["status"], "failed")
            self.assertEqual(writer.metadata["frame_count"], 0)

    def test_source_change_is_detected_and_partial_data_is_not_completed(self):
        with tempfile.TemporaryDirectory() as folder:
            source = Path(folder)/"sample.bonds"
            source.write_text(_frame(0, {}, 1))
            changed = []
            def progress(_):
                if not changed:
                    with source.open("a") as handle:
                        handle.write("\n# source changed\n")
                    changed.append(True)
            with self.assertRaisesRegex(ValueError, "Source file changed"):
                analyze_counts_file(source, ELEMENTS, Path(folder)/"changed", progress_callback=progress)
            metadata = json.loads((Path(folder)/"changed/count_analysis.json").read_text())
            self.assertEqual(metadata["status"], "failed")


if __name__ == "__main__":
    unittest.main()
