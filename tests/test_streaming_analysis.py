from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from reaxnet.analysis import AnalyzerConfig
from reaxnet.streaming import analyze_file_streaming


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _three_atom_frame(timestep: int, edge: tuple[int, int]) -> str:
    adjacency = {1: [], 2: [], 3: []}
    left, right = edge
    adjacency[left].append(right)
    adjacency[right].append(left)
    lines = [
        f"# Timestep {timestep}",
        "# Number of particles 3",
        "# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q",
    ]
    for atom_id in (1, 2, 3):
        atom_type = 1 if atom_id == 1 else 2
        neighbors = adjacency[atom_id]
        if neighbors:
            lines.append(
                f"{atom_id} {atom_type} 1 {neighbors[0]} 0 0.950 0.950 0.000 0.000"
            )
        else:
            lines.append(f"{atom_id} {atom_type} 0 0 0.000 0.000 0.000")
    return "\n".join(lines)


class StreamingAnalysisTests(unittest.TestCase):
    def test_streaming_outputs_match_small_trajectory_contract(self) -> None:
        progress: list[dict[str, int | float]] = []
        with tempfile.TemporaryDirectory() as tmp:
            result = analyze_file_streaming(
                PROJECT_ROOT / "examples" / "ch4_dissociation.bonds",
                {1: "C", 2: "H"},
                tmp,
                AnalyzerConfig(timestep_scale=0.001),
                progress_callback=progress.append,
            )

            self.assertEqual(result.summary["frame_count"], 2)
            self.assertEqual(result.summary["species_count"], 3)
            self.assertEqual(result.summary["accepted_event_count"], 1)
            self.assertTrue(result.summary["streaming_analysis"])
            self.assertEqual(result.network["reactions"][0]["net_count"], 1)
            self.assertEqual(progress[-1]["bytes_processed"], progress[-1]["bytes_total"])
            with result.paths["reaction_events"].open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["reaction"], "C -> [CH3] + [H]")

    def test_multistep_recrossing_is_preserved_in_streaming_mode(self) -> None:
        trajectory = "\n".join(
            [
                _three_atom_frame(0, (1, 2)),
                _three_atom_frame(10, (1, 3)),
                _three_atom_frame(20, (2, 3)),
                _three_atom_frame(30, (1, 2)),
            ]
        )
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "cycle.bonds"
            source.write_text(trajectory, encoding="utf-8")
            result = analyze_file_streaming(
                source,
                {1: "C", 2: "H"},
                Path(tmp) / "out",
                AnalyzerConfig(
                    timestep_scale=0.001,
                    recrossing_window=0.04,
                    generate_identifiers=False,
                ),
            )

            self.assertEqual(result.summary["detected_event_count"], 3)
            self.assertEqual(result.summary["filtered_event_count"], 3)
            self.assertEqual(result.summary["accepted_event_count"], 0)
            self.assertEqual(result.network["reactions"], [])
            with result.paths["reaction_events"].open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual({row["filter_group"] for row in rows}, {"RC000001"})
            self.assertEqual(
                {row["filter_reason"] for row in rows},
                {"multistep_state_return_within_window"},
            )

    def test_static_long_trajectory_keeps_one_state_in_memory(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "static.bonds"
            with source.open("w", encoding="utf-8") as handle:
                for timestep in range(2000):
                    handle.write(
                        f"# Timestep {timestep}\n"
                        "# Number of particles 1\n"
                        "1 1 0 0 0.000 0.000 0.000\n"
                    )
            result = analyze_file_streaming(
                source,
                {1: "H"},
                Path(tmp) / "out",
                AnalyzerConfig(
                    timestep_scale=0.001,
                    recrossing_window=0.1,
                    generate_identifiers=False,
                ),
            )

            self.assertEqual(result.summary["frame_count"], 2000)
            self.assertEqual(result.summary["detected_event_count"], 0)
            self.assertEqual(result.summary["peak_buffered_states"], 1)
            self.assertEqual(result.summary["peak_pending_events"], 0)


if __name__ == "__main__":
    unittest.main()
