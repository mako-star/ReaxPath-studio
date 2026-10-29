"""Bounded CHO bond-type and electron-state search on fixed observed graphs."""
from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass

from reaxnet.identity import CandidateIdentifiers, identify_structure_candidate

MODEL_VERSION = 'cho-local-states-1'
NORMALIZATION_VERSION = 'rdkit-aromatic-smiles-1-no-tautomer-or-charge-normalization'

# (sum of integer bond orders, formal charge, number of radical electrons).
# This is an explicitly limited model, not an exhaustive electronic-state model.
# In particular it excludes singlet carbenes and does not resolve excited states or spin multiplicities.
LOCAL_STATES = {
    'H': ((1, 0, 0), (0, 0, 1), (0, 1, 0), (0, -1, 0)),
    'C': ((4, 0, 0), (3, 0, 1), (2, 0, 2), (3, 1, 0), (3, -1, 0)),
    'O': ((2, 0, 0), (1, 0, 1), (0, 0, 2), (1, -1, 0), (3, 1, 0)),
}


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def rdkit_version():
    try:
        from rdkit import rdBase
        return rdBase.rdkitVersion
    except ImportError:
        return 'unavailable'


@dataclass(frozen=True)
class StructureConfig:
    total_charges: tuple[int, ...]
    max_radicals: int
    allow_charge_separation: bool = False
    max_search_nodes: int = 50000
    max_candidates: int = 128
    max_atoms: int = 128
    max_variable_edges: int = 128
    cache_size: int = 256

    def __post_init__(self):
        if not self.total_charges or any(type(q) is not int for q in self.total_charges):
            raise ValueError('total_charges must explicitly contain integer fragment charges')
        object.__setattr__(self, 'total_charges', tuple(sorted(set(self.total_charges))))
        if type(self.allow_charge_separation) is not bool:
            raise ValueError('allow_charge_separation must be boolean')
        if type(self.max_radicals) is not int or self.max_radicals < 0:
            raise ValueError('max_radicals must be a nonnegative integer')
        for name in ('max_search_nodes', 'max_candidates', 'max_atoms', 'max_variable_edges', 'cache_size'):
            value = getattr(self, name)
            if type(value) is not int or value < (0 if name == 'cache_size' else 1):
                raise ValueError(f'{name} has an invalid limit')
        if self.max_atoms > 256 or self.max_variable_edges > 256:
            raise ValueError('recursive search supports at most 256 atoms and 256 variable edges')

    @property
    def assumption_id(self):
        return 'A' + digest({'config': asdict(self), 'model': MODEL_VERSION,
                            'normalization': NORMALIZATION_VERSION, 'rdkit': rdkit_version()})


@dataclass(frozen=True)
class MolecularGraph:
    symbols: tuple[str, ...]
    edges: tuple[tuple[int, int], ...]

    def __post_init__(self):
        object.__setattr__(self, 'symbols', tuple(self.symbols))
        if not self.symbols:
            raise ValueError('empty molecule graph')
        if any(len(e) != 2 or any(type(a) is not int for a in e) for e in self.edges):
            raise ValueError('edges require pairs of integer local indices')
        edges = tuple(sorted(tuple(sorted(e)) for e in self.edges))
        if len(set(edges)) != len(edges) or any(len(e) != 2 or e[0] < 0 or e[1] >= len(self.symbols) or e[0] == e[1] for e in edges):
            raise ValueError('invalid molecule edges')
        object.__setattr__(self, 'edges', edges)
        seen, pending = {0}, [0]
        adjacency = [[] for _ in self.symbols]
        for a, b in edges:
            adjacency[a].append(b)
            adjacency[b].append(a)
        while pending:
            for b in adjacency[pending.pop()]:
                if b not in seen:
                    seen.add(b)
                    pending.append(b)
        if len(seen) != len(self.symbols):
            raise ValueError('a molecule graph must be connected')

    @classmethod
    def from_instance(cls, instance, symbols_by_atom):
        positions = {a: i for i, a in enumerate(instance.atom_ids)}
        return cls(tuple(symbols_by_atom[a] for a in instance.atom_ids),
                   tuple((positions[a], positions[b]) for a, b in instance.connections))


@dataclass(frozen=True)
class StructureCandidate:
    bond_orders: tuple[int, ...]
    formal_charges: tuple[int, ...]
    radical_electrons: tuple[int, ...]
    identifiers: CandidateIdentifiers


@dataclass(frozen=True)
class InferenceResult:
    identity_status: str
    search_complete: bool
    reason: str
    assumption_id: str
    search_nodes: int
    candidates: tuple[StructureCandidate, ...]


class _SearchLimit(Exception):
    pass


def infer_structure(graph: MolecularGraph, config: StructureConfig) -> InferenceResult:
    assumption = config.assumption_id
    candidates = []
    nodes = 0
    def result(status, complete, reason):
        return InferenceResult(status, complete, reason, assumption, nodes, tuple(candidates))
    if set(graph.symbols) - LOCAL_STATES.keys():
        return result('unresolved', False, 'unsupported_elements')
    if len(graph.symbols) > config.max_atoms:
        return result('unresolved', False, 'max_atoms_exceeded')
    domains = []
    for symbol in graph.symbols:
        states = [s for s in LOCAL_STATES[symbol] if s[2] <= config.max_radicals]
        if not config.allow_charge_separation:
            states = [s for s in states if s[1] == 0 or any(q * s[1] > 0 for q in config.total_charges)]
        domains.append(states)
    valences = [0] * len(graph.symbols)
    orders = [1] * len(graph.edges)
    variable = []
    remaining = [0] * len(graph.symbols)
    for index, (a, b) in enumerate(graph.edges):
        valences[a] += 1
        valences[b] += 1
        if graph.symbols[a] != 'H' and graph.symbols[b] != 'H':
            variable.append(index)
            remaining[a] += 2
            remaining[b] += 2
    if len(variable) > config.max_variable_edges:
        return result('unresolved', False, 'max_variable_edges_exceeded')

    def tick():
        nonlocal nodes
        if nodes >= config.max_search_nodes:
            raise _SearchLimit('max_search_nodes_exceeded')
        nodes += 1

    def feasible():
        compatible = [[s for s in states if valences[i] <= s[0] <= valences[i] + remaining[i]]
                      for i, states in enumerate(domains)]
        if any(not states for states in compatible):
            return False
        if sum(min(s[2] for s in states) for states in compatible) > config.max_radicals:
            return False
        low = sum(min(s[1] for s in states) for states in compatible)
        high = sum(max(s[1] for s in states) for states in compatible)
        return any(low <= q <= high for q in config.total_charges)

    def assign_states():
        choices = [[s for s in states if s[0] == valences[i]] for i, states in enumerate(domains)]
        charge_min = [0] * (len(choices) + 1)
        charge_max = [0] * (len(choices) + 1)
        radical_min = [0] * (len(choices) + 1)
        for i in reversed(range(len(choices))):
            charge_min[i] = charge_min[i + 1] + min(s[1] for s in choices[i])
            charge_max[i] = charge_max[i + 1] + max(s[1] for s in choices[i])
            radical_min[i] = radical_min[i + 1] + min(s[2] for s in choices[i])
        charges, radicals = [], []
        def visit(i, q, r, absolute_q):
            tick()
            if r + radical_min[i] > config.max_radicals:
                return
            if not any(q + charge_min[i] <= target <= q + charge_max[i] for target in config.total_charges):
                return
            if i == len(choices):
                if not config.allow_charge_separation and absolute_q != abs(q):
                    return
                if len(candidates) >= config.max_candidates:
                    raise _SearchLimit('max_candidates_exceeded')
                identity = identify_structure_candidate(graph.symbols, graph.edges, tuple(orders), tuple(charges), tuple(radicals))
                candidates.append(StructureCandidate(tuple(orders), tuple(charges), tuple(radicals), identity))
                return
            for _, charge, radical in choices[i]:
                charges.append(charge)
                radicals.append(radical)
                visit(i + 1, q + charge, r + radical, absolute_q + abs(charge))
                charges.pop()
                radicals.pop()
        visit(0, 0, 0, 0)

    def assign_bonds(i):
        tick()
        if not feasible():
            return
        if i == len(variable):
            assign_states()
            return
        index = variable[i]
        a, b = graph.edges[index]
        remaining[a] -= 2
        remaining[b] -= 2
        for extra in range(3):
            valences[a] += extra
            valences[b] += extra
            orders[index] = 1 + extra
            assign_bonds(i + 1)
            valences[a] -= extra
            valences[b] -= extra
        remaining[a] += 2
        remaining[b] += 2
        orders[index] = 1

    try:
        assign_bonds(0)
    except _SearchLimit as exc:
        return result('unresolved', False, str(exc))
    if not candidates:
        return result('unresolved', True, 'no_candidate_under_declared_model')
    if any(c.identifiers.smiles_status != 'ok' for c in candidates):
        return result('unresolved', True, 'candidate_validation_or_smiles_incomplete')
    unique = {c.identifiers.canonical_smiles for c in candidates}
    if len(unique) == 1:
        keys = {c.identifiers.inchikey for c in candidates if c.identifiers.inchi_status == 'ok'}
        if len(keys) > 1:
            return result('unresolved', True, 'identifier_disagreement_for_same_smiles')
        return result('supported', True, 'unique_standardized_structure_under_declared_model')
    return result('ambiguous', True, 'multiple_standardized_structures')
