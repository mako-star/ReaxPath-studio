from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from reaxnet.analysis import AnalyzerConfig, analyze_frames
from reaxnet.chemistry import build_frame_state
from reaxnet.outputs import write_outputs
from reaxnet.parser import parse_bonds_lines, parse_element_mapping


ELEMENTS = {1: "C", 2: "H"}


def _frame(
    timestep: int,
    edges: dict[tuple[int, int], float],
    atom_count: int = 5,
    atom_types: dict[int, int] | None = None,
) -> str:
    normalized = {(min(a, b), max(a, b)): order for (a, b), order in edges.items()}
    adjacency: dict[int, list[tuple[int, float]]] = {atom_id: [] for atom_id in range(1, atom_count + 1)}
    for (a, b), order in normalized.items():
        adjacency[a].append((b, order))
        adjacency[b].append((a, order))

    lines = [
        f"# Timestep {timestep}",
        "#",
        f"# Number of particles {atom_count}",
        "#",
        "# Max number of bonds per atom 8 with coarse bond order cutoff 0.300",
        "# Particle connection table and bond orders",
        "# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q",
    ]
    for atom_id in range(1, atom_count + 1):
        neighbors = sorted(adjacency[atom_id])
        atom_type = (atom_types or {}).get(atom_id, 1 if atom_id == 1 else 2)
        parts = [str(atom_id), str(atom_type), str(len(neighbors))]
        parts.extend(str(neighbor) for neighbor, _ in neighbors)
        parts.append("0")
        parts.extend(f"{order:.3f}" for _, order in neighbors)
        parts.extend((f"{sum(order for _, order in neighbors):.3f}", "0.000", "0.000"))
        lines.append(" ".join(parts))
    lines.append("#")
    return "\n".join(lines)


def _trajectory(*frames: str) -> str:
    return "\n".join(frames) + "\n"


CH4_BONDS = {(1, 2): 0.96, (1, 3): 0.95, (1, 4): 0.94, (1, 5): 0.93}
CH3_H_BONDS = {(1, 2): 0.96, (1, 3): 0.95, (1, 4): 0.94}


class ParserTests(unittest.TestCase):
    def test_parses_standard_fix_reax_bonds_blocks(self) -> None:
        text = _trajectory(_frame(0, CH4_BONDS), _frame(10, CH3_H_BONDS))

        frames = list(parse_bonds_lines(text.splitlines()))

        self.assertEqual([frame.timestep for frame in frames], [0, 10])
        self.assertEqual(frames[0].expected_atom_count, 5)
        self.assertEqual(frames[0].atoms[1].neighbors, (2, 3, 4, 5))
        self.assertEqual(frames[0].atoms[1].bond_orders, (0.96, 0.95, 0.94, 0.93))
        self.assertAlmostEqual(frames[0].atoms[1].charge, 0.0)

    def test_parses_type_to_element_mapping(self) -> None:
        self.assertEqual(parse_element_mapping("1:C, 2:H,3:O"), {1: "C", 2: "H", 3: "O"})


class ChemistryTests(unittest.TestCase):
    def test_file_topology_keeps_fifth_carbon_bond(self) -> None:
        five_bonds = {
            (1, 2): 1.00,
            (1, 3): 0.99,
            (1, 4): 0.98,
            (1, 5): 0.97,
            (1, 6): 0.96,
        }
        frame = next(parse_bonds_lines(_frame(0, five_bonds, atom_count=6).splitlines()))

        state = build_frame_state(frame, ELEMENTS, generate_identifiers=False)

        self.assertEqual(len(state.bonds), 5)
        self.assertEqual(len(state.transition_bonds), 0)
        self.assertEqual(sorted(molecule.formula for molecule in state.molecules), ["CH5"])


class AnalysisTests(unittest.TestCase):
    def test_bo_crossing_old_thresholds_does_not_create_events(self) -> None:
        frames = parse_bonds_lines(_trajectory(*(
            _frame(i * 10, {(1, 2): order}, atom_count=2)
            for i, order in enumerate((0.9, 0.5, 0.35, 0.3, 0.0))
        ), _frame(50, {}, atom_count=2)).splitlines())
        result = analyze_frames(frames, ELEMENTS, AnalyzerConfig(generate_identifiers=False))
        self.assertEqual(len(result.events), 1)
        self.assertEqual(result.events[0].timestep, 50)
        self.assertEqual(result.events[0].broken_bonds, ((1, 2),))
        self.assertEqual(result.network['metadata']['connectivity_policy'], 'file_neighbor_list')

    def test_detects_dissociation_and_generates_identifiers(self) -> None:
        frames = parse_bonds_lines(
            _trajectory(_frame(0, CH4_BONDS), _frame(10, CH3_H_BONDS)).splitlines()
        )

        result = analyze_frames(
            frames,
            ELEMENTS,
            AnalyzerConfig(timestep_scale=0.001),
        )

        self.assertEqual(len(result.events), 1)
        event = result.events[0]
        self.assertTrue(event.accepted)
        self.assertEqual(event.event_type, "dissociation")
        self.assertEqual(event.reactant_formulas, ("CH4",))
        self.assertEqual(sorted(event.product_formulas), ["CH3", "H"])
        self.assertEqual(event.broken_bonds, ((1, 5),))
        self.assertEqual(event.formed_bonds, ())

        methane = next(entry for entry in result.registry.values() if entry.formula == "CH4")
        self.assertEqual(methane.smiles, "C")
        self.assertTrue(methane.inchi.startswith("InChI=1S/CH4"))
        self.assertEqual(len(methane.inchikey), 27)
        self.assertEqual(len(result.network["reactions"]), 1)

    def test_exact_reverse_inside_window_is_filtered_as_recrossing(self) -> None:
        frames = parse_bonds_lines(
            _trajectory(
                _frame(0, CH4_BONDS),
                _frame(10, CH3_H_BONDS),
                _frame(20, CH4_BONDS),
            ).splitlines()
        )

        result = analyze_frames(
            frames,
            ELEMENTS,
            AnalyzerConfig(
                timestep_scale=0.001,
                recrossing_window=0.02,
            ),
        )

        self.assertEqual(len(result.events), 2)
        self.assertTrue(all(not event.accepted for event in result.events))
        self.assertTrue(all(event.filter_reason == "exact_reverse_within_window" for event in result.events))
        self.assertEqual({event.filter_group for event in result.events}, {"RC000001"})
        self.assertEqual(result.network["reactions"], [])

    def test_multistep_atom_path_return_filters_entire_event_chain(self) -> None:
        atom_types = {1: 1, 2: 2, 3: 2}
        frames = parse_bonds_lines(
            _trajectory(
                _frame(0, {(1, 2): 0.95}, atom_count=3, atom_types=atom_types),
                _frame(10, {(1, 3): 0.95}, atom_count=3, atom_types=atom_types),
                _frame(20, {(2, 3): 0.95}, atom_count=3, atom_types=atom_types),
                _frame(30, {(1, 2): 0.95}, atom_count=3, atom_types=atom_types),
            ).splitlines()
        )

        result = analyze_frames(
            frames,
            ELEMENTS,
            AnalyzerConfig(
                timestep_scale=0.001,
                recrossing_window=0.04,
                generate_identifiers=False,
            ),
        )

        self.assertEqual(len(result.events), 3)
        self.assertTrue(all(not event.accepted for event in result.events))
        self.assertEqual(
            {event.filter_reason for event in result.events},
            {"multistep_state_return_within_window"},
        )
        self.assertEqual({event.filter_group for event in result.events}, {"RC000001"})
        self.assertEqual(result.network["reactions"], [])

    def test_writes_complete_output_contract(self) -> None:
        frames = parse_bonds_lines(
            _trajectory(_frame(0, CH4_BONDS), _frame(10, CH3_H_BONDS)).splitlines()
        )
        result = analyze_frames(frames, ELEMENTS, AnalyzerConfig(timestep_scale=0.001))

        with tempfile.TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            written = write_outputs(result, output_dir)

            self.assertEqual(
                set(written),
                {
                    "species_timeseries",
                    "reaction_events",
                    "molecule_registry",
                    "reaction_network",
                },
            )
            self.assertTrue(all(path.is_file() for path in written.values()))

            with written["molecule_registry"].open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual({row["formula"] for row in rows}, {"CH4", "CH3", "H"})
            self.assertTrue(all("status" in row and "message" in row for row in rows))


if __name__ == "__main__":
    unittest.main()
