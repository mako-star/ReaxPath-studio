from __future__ import annotations

import unittest
from unittest.mock import patch

from reaxnet.chemistry import TopologyRegistry, build_frame_state
from reaxnet.connectivity import ConnectivityConfig, ConnectivityTracker
from reaxnet.identity import IdentityResult
from reaxnet.models import AtomRecord, Bond, Frame, FrameState, MoleculeRecord
from reaxnet.parser import BondsParseError, parse_bonds_lines


ELEMENTS = {1: "H", 2: "C", 3: "O"}
EDGES = {(3, 80): 0.41, (7, 12): 1.2, (12, 40): 0.76}


def sparse_frame(timestep=40, edges=None):
    # Input order deliberately differs from both atom ID and molecule order.
    types = {80: 3, 7: 1, 12: 3, 3: 1, 99: 2, 40: 1}
    neighbors = {atom_id: [] for atom_id in types}
    for (a, b), order in (EDGES if edges is None else edges).items():
        neighbors[a].append((b, order))
        neighbors[b].append((a, order))
    atoms = {
        atom_id: AtomRecord(
            atom_id, atom_type,
            tuple(n for n, _ in neighbors[atom_id]), 0,
            tuple(order for _, order in neighbors[atom_id]),
            sum(order for _, order in neighbors[atom_id]), 0.0, 0.0,
        )
        for atom_id, atom_type in types.items()
    }
    return Frame(timestep, atoms, expected_atom_count=len(atoms))


class FrameStateOptimizationTests(unittest.TestCase):
    def test_sparse_unsorted_ids_match_preoptimization_state(self):
        # Golden state captured before the optimization. In particular, the
        # tracker F... IDs must not replace the existing public T... IDs.
        disabled = {"status": "disabled", "message": "Identifier generation was disabled"}
        expected = FrameState(
            timestep=40,
            bonds={
                (3, 80): Bond(3, 80, 0.41, 0.5),
                (7, 12): Bond(7, 12, 1.2, 1.0),
                (12, 40): Bond(12, 40, 0.76, 1.0),
            },
            transition_bonds=(),
            molecules=(
                MoleculeRecord((3, 80), "HO", "HO-b95c0829f7ea16ff", "1-2:0.5", 2, 1,
                    instance_id="T40:M00000000", structure_id="HO-001",
                    atom_elements=("H", "O"), connections=((3, 80),), **disabled),
                MoleculeRecord((7, 12, 40), "H2O", "H2O-330bafe7c7cb5b84", "1-2:1;2-3:1", 3, 2,
                    instance_id="T40:M00000001", structure_id="H2O-001",
                    atom_elements=("H", "O", "H"), connections=((7, 12), (12, 40)), **disabled),
                MoleculeRecord((99,), "C", "C-d86ee1257bb45818", "", 1, 0,
                    instance_id="T40:M00000002", structure_id="C-001",
                    atom_elements=("C",), connections=(), **disabled),
            ),
            atom_to_molecule={3: 0, 80: 0, 7: 1, 12: 1, 40: 1, 99: 2},
            symbols={80: "O", 7: "H", 12: "O", 3: "H", 99: "C", 40: "H"},
        )
        frame = sparse_frame()
        observed = ConnectivityTracker(ELEMENTS, ConnectivityConfig()).process(frame)
        for supplied in (None, observed):
            with self.subTest(precomputed_connectivity=supplied is not None):
                actual = build_frame_state(frame, ELEMENTS, generate_identifiers=False,
                                           connectivity_frame=supplied)
                self.assertEqual(actual, expected)

    def test_removed_edges_are_not_restored_from_tracker_history(self):
        tracker = ConnectivityTracker(ELEMENTS, ConnectivityConfig(expected_interval=10))
        registry = TopologyRegistry()
        first = sparse_frame()
        before = build_frame_state(first, ELEMENTS, generate_identifiers=False,
                                   connectivity_frame=tracker.process(first), topology_registry=registry)
        second = sparse_frame(50, {(7, 12): 1.2, (12, 40): 0.76})
        after = build_frame_state(second, ELEMENTS, generate_identifiers=False,
                                  connectivity_frame=tracker.process(second), topology_registry=registry)
        self.assertEqual(set(after.bonds), {(7, 12), (12, 40)})
        self.assertEqual([m.atom_ids for m in after.molecules], [(3,), (7, 12, 40), (80,), (99,)])
        self.assertEqual([m.instance_id for m in after.molecules],
                         [f"T50:M{i:08d}" for i in range(4)])
        self.assertEqual(after.molecules[1].species_key, before.molecules[1].species_key)
        self.assertEqual(after.molecules[1].structure_id, "H2O-001")

    def test_disabled_identifiers_still_preserve_changing_bond_signature(self):
        registry = TopologyRegistry()
        states = [build_frame_state(sparse_frame(i * 10, {(3, 80): order}), ELEMENTS,
                                    generate_identifiers=False, topology_registry=registry)
                  for i, order in enumerate((0.0, 1.3, 2.0))]
        self.assertEqual([s.molecules[0].bond_signature for s in states],
                         ["1-2:0", "1-2:1.5", "1-2:2"])
        self.assertEqual(len({s.molecules[0].species_key for s in states}), 1)
        for state in states:
            molecule = state.molecules[0]
            self.assertEqual(molecule.connections, ((3, 80),))
            self.assertEqual((molecule.status, molecule.message),
                             ("disabled", "Identifier generation was disabled"))
            self.assertEqual((molecule.smiles, molecule.inchi, molecule.inchikey), ("", "", ""))

    def test_enabled_cache_distinguishes_bond_interpretations_and_can_reuse_them(self):
        cache = {}
        registry = TopologyRegistry()
        # Other atoms are isolated; give them deterministic identities as well.
        def identify(atom_ids, symbols, bonds, *, enabled):
            if not enabled:
                return IdentityResult(status="disabled", message="Identifier generation was disabled")
            value = bonds.get((3, 80))
            return IdentityResult(smiles=f"order-{value.effective_order}" if value else "isolated", status="ok")

        with patch("reaxnet.chemistry.identify_molecule", side_effect=identify):
            states = [build_frame_state(sparse_frame(i * 10, {(3, 80): order}), ELEMENTS,
                                        identity_cache=cache, topology_registry=registry)
                      for i, order in enumerate((0.4, 2.0, 0.4))]
            disabled = build_frame_state(sparse_frame(30, {(3, 80): 0.4}), ELEMENTS,
                generate_identifiers=False, identity_cache=cache, topology_registry=registry)
        self.assertEqual([s.molecules[0].smiles for s in states], ["order-0.5", "order-2.0", "order-0.5"])
        self.assertEqual(len({s.molecules[0].structure_id for s in states}), 1)
        self.assertEqual(disabled.molecules[0].status, "disabled")
        self.assertEqual(disabled.molecules[0].smiles, "")
        self.assertEqual(disabled.molecules[0].bond_signature, states[0].molecules[0].bond_signature)


class ParserOptimizationTests(unittest.TestCase):
    def test_missing_neighbors_report_first_input_atom_and_complete_sorted_set(self):
        text = "\n".join([
            "# Timestep 40", "# Number of particles 2",
            "99 1 2 888 777 0 0.7 0.8 1.5 0 0",
            "2 1 1 666 0 0.5 0.5 0 0",
        ])
        for suffix in ("", "\n# Timestep 50\n# Number of particles 0"):
            with self.subTest(next_frame=bool(suffix)), self.assertRaises(BondsParseError) as caught:
                list(parse_bonds_lines((text + suffix).splitlines()))
            self.assertEqual(str(caught.exception),
                "Timestep 40: atom 99 references missing atoms [777, 888]")

    def test_explicit_empty_frame_remains_valid_and_undeclared_empty_frame_fails(self):
        frame = next(parse_bonds_lines(["# Timestep 0", "# Number of particles 0"]))
        observed = ConnectivityTracker(ELEMENTS, ConnectivityConfig()).process(frame)
        for supplied in (None, observed):
            with self.subTest(precomputed_connectivity=supplied is not None):
                state = build_frame_state(frame, ELEMENTS, generate_identifiers=False,
                                          connectivity_frame=supplied)
                self.assertEqual(state, FrameState(0, {}, (), (), {}, {}))
        with self.assertRaisesRegex(BondsParseError, "empty frame without an explicit zero atom count"):
            list(parse_bonds_lines(["# Timestep 0"]))


if __name__ == "__main__":
    unittest.main()
