from __future__ import annotations

import hashlib
import math
from collections import Counter, defaultdict, deque
from collections.abc import Mapping, MutableMapping

from reaxnet.identity import IdentityResult, identify_molecule
from reaxnet.models import Bond, BondKey, Frame, FrameState, MoleculeRecord


_DISABLED_IDENTITY = IdentityResult(
    status="disabled", message="Identifier generation was disabled"
)


class ChemistryError(ValueError):
    """Raised when atom typing is insufficient for analysis."""


def _half_step_bond_order(raw_order: float) -> float:
    return math.floor(raw_order * 2.0 + 0.5) / 2.0


def _symbols_for_frame(frame: Frame, element_mapping: Mapping[int, str]) -> dict[int, str]:
    symbols: dict[int, str] = {}
    for atom_id, atom in frame.atoms.items():
        try:
            symbols[atom_id] = element_mapping[atom.atom_type]
        except KeyError as exc:
            raise ChemistryError(
                f"Timestep {frame.timestep}: LAMMPS atom type {atom.atom_type} "
                f"(atom {atom_id}) has no --elements mapping"
            ) from exc
    return symbols


def select_real_bonds(
    frame: Frame,
    element_mapping: Mapping[int, str],
    observed=None,
) -> tuple[dict[BondKey, Bond], tuple[Bond, ...]]:
    """Use validated file adjacency without BO or valence-based edge removal.

    effective_order is retained only for the legacy identity representation.
    The PRD structure pipeline infers chemical bond types separately.
    """
    from reaxnet.connectivity import ConnectivityConfig, ConnectivityTracker

    if observed is None:
        observed = ConnectivityTracker(element_mapping, ConnectivityConfig()).process(frame)
    bonds = {}
    for pair in observed.observations:
        if pair.connected:
            order = pair.merged_order
            bonds[pair.atom_ids] = Bond(pair.atom_ids[0], pair.atom_ids[1], order,
                                        _half_step_bond_order(order))
    return bonds, ()


def _connected_components(atom_ids: set[int], bonds: Mapping[BondKey, Bond]) -> list[tuple[int, ...]]:
    adjacency: dict[int, set[int]] = {atom_id: set() for atom_id in atom_ids}
    for atom1, atom2 in bonds:
        adjacency[atom1].add(atom2)
        adjacency[atom2].add(atom1)

    components: list[tuple[int, ...]] = []
    unseen = set(atom_ids)
    while unseen:
        start = min(unseen)
        queue = deque([start])
        unseen.remove(start)
        component: list[int] = []
        while queue:
            atom_id = queue.popleft()
            component.append(atom_id)
            for neighbor in sorted(adjacency[atom_id]):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    queue.append(neighbor)
        components.append(tuple(sorted(component)))
    return components


def molecular_formula(atom_ids: tuple[int, ...], symbols: Mapping[int, str]) -> str:
    counts = Counter(symbols[atom_id] for atom_id in atom_ids)
    if "C" in counts:
        order = ["C"]
        if "H" in counts:
            order.append("H")
        order.extend(sorted(symbol for symbol in counts if symbol not in {"C", "H"}))
    else:
        order = sorted(counts)
    return "".join(
        symbol + (str(counts[symbol]) if counts[symbol] != 1 else "")
        for symbol in order
    )


def _graph_fingerprint(
    atom_ids: tuple[int, ...],
    symbols: Mapping[int, str],
    bonds: Mapping[BondKey, Bond],
) -> str:
    adjacency: dict[int, list[int]] = {atom_id: [] for atom_id in atom_ids}
    for atom1, atom2 in bonds:
        if atom1 in adjacency and atom2 in adjacency:
            adjacency[atom1].append(atom2)
            adjacency[atom2].append(atom1)

    labels = {
        atom_id: hashlib.sha256(symbols[atom_id].encode("ascii")).hexdigest()
        for atom_id in atom_ids
    }
    for _ in range(max(2, min(len(atom_ids), 8))):
        labels = {
            atom_id: hashlib.sha256(
                (labels[atom_id] + "|" + "|".join(sorted(labels[n] for n in adjacency[atom_id]))).encode(
                    "ascii"
                )
            ).hexdigest()
            for atom_id in atom_ids
        }
    payload = "|".join(sorted(labels.values())) + f"|edges={sum(map(len, adjacency.values())) // 2}"
    return hashlib.sha256(payload.encode("ascii")).hexdigest()[:16]


def _bond_signature(atom_ids: tuple[int, ...], bonds: Mapping[BondKey, Bond]) -> str:
    local_index = {atom_id: index for index, atom_id in enumerate(atom_ids, start=1)}
    return ";".join(
        f"{local_index[atom1]}-{local_index[atom2]}:{bond.effective_order:g}"
        for (atom1, atom2), bond in sorted(bonds.items())
        if atom1 in local_index and atom2 in local_index
    )


class TopologyRegistry:
    """Run-local exact topology registry; fingerprints only choose comparison buckets."""
    def __init__(self):
        self.buckets = defaultdict(list)
        self.exact = {}
        self.structure_ids = {}
        self.formula_counts = Counter()

    def structure_id(self, species_key, formula):
        if species_key not in self.structure_ids:
            self.formula_counts[formula] += 1
            self.structure_ids[species_key] = f"{formula}-{self.formula_counts[formula]:03d}"
        return self.structure_ids[species_key]

    def register(self, component, symbols, bonds, formula):
        indices = {a:i for i,a in enumerate(component)}
        labels = tuple(symbols[a] for a in component)
        edges = tuple(sorted((indices[a],indices[b]) for a,b in bonds))
        signature = (labels,edges)
        if signature in self.exact:
            return self.exact[signature]
        bucket_key = f"{formula}-{_graph_fingerprint(component,symbols,bonds)}"
        bucket = self.buckets[bucket_key]
        for old_signature,key in bucket:
            if _isomorphic(signature,old_signature):
                self.exact[signature] = key
                return key
        key = bucket_key if not bucket else f"{bucket_key}-{len(bucket)+1}"
        bucket.append((signature,key))
        self.exact[signature] = key
        return key


def _isomorphic(left, right, max_nodes=200000):
    """Exact adjacency mapping. Exhausted budget conservatively avoids merging."""
    ls,le=left;rs,re=right
    if len(ls)!=len(rs) or len(le)!=len(re): return False
    la=[set() for _ in ls];ra=[set() for _ in rs]
    for a,b in le: la[a].add(b);la[b].add(a)
    for a,b in re: ra[a].add(b);ra[b].add(a)
    signatures=lambda sy,adj:[(sy[i],len(adj[i]),tuple(sorted(sy[j] for j in adj[i]))) for i in range(len(sy))]
    lp,rp=signatures(ls,la),signatures(rs,ra)
    if Counter(lp)!=Counter(rp): return False
    by_signature = defaultdict(list)
    for b, signature in enumerate(rp):
        by_signature[signature].append(b)
    choices={a:by_signature[lp[a]] for a in range(len(ls))}
    mapping={};used=set();stack=[];nodes=0
    while len(mapping)<len(ls):
        if len(stack)==len(mapping):
            a=min((a for a in choices if a not in mapping),key=lambda a:(-len(la[a]&mapping.keys()),len(choices[a]),-len(la[a]),a))
            mapped_neighbors = [mapping[v] for v in la[a] if v in mapping]
            if mapped_neighbors:
                possible = set(ra[mapped_neighbors[0]])
                for neighbor in mapped_neighbors[1:]:
                    possible.intersection_update(ra[neighbor])
                candidates = sorted(b for b in possible if rp[b] == lp[a] and b not in used)
            else:
                candidates = choices[a]
            stack.append((a,iter(candidates)))
        a,candidates=stack[-1]
        for b in candidates:
            nodes+=1
            if nodes>max_nodes: return False
            if b not in used and all((v in la[a])==(w in ra[b]) for v,w in mapping.items()):
                mapping[a]=b;used.add(b);break
        else:
            stack.pop()
            if not stack: return False
            previous=stack[-1][0];used.remove(mapping.pop(previous))
    return True


def build_frame_state(
    frame: Frame,
    element_mapping: Mapping[int, str],
    *,
    generate_identifiers: bool = True,
    identity_cache: MutableMapping[object, IdentityResult] | None = None,
    connectivity_frame=None,
    topology_registry: TopologyRegistry | None = None,
) -> FrameState:
    symbols = _symbols_for_frame(frame, element_mapping)
    if connectivity_frame is None:
        from reaxnet.connectivity import ConnectivityConfig, ConnectivityTracker

        connectivity_frame = ConnectivityTracker(element_mapping, ConnectivityConfig()).process(frame)
    bonds, transition_bonds = select_real_bonds(
        frame,
        element_mapping,
        connectivity_frame,
    )
    topology_registry = topology_registry or TopologyRegistry()
    molecules: list[MoleculeRecord] = []
    atom_to_molecule: dict[int, int] = {}

    # The tracker already partitions atoms and edges in the same stable order.
    # Reuse that work instead of recomputing components and scanning all frame
    # edges separately for every molecule.
    for instance in connectivity_frame.instances:
        component = instance.atom_ids
        component_bonds = {key: bonds[key] for key in instance.connections}
        formula = molecular_formula(component, symbols)
        species_key = topology_registry.register(component, symbols, component_bonds, formula)
        atom_elements = tuple(symbols[a] for a in component)
        bond_signature = _bond_signature(component, component_bonds)
        if generate_identifiers:
            identity_key = (atom_elements, bond_signature, generate_identifiers)
            identity = identity_cache.get(identity_key) if identity_cache is not None else None
            if identity is None:
                identity = identify_molecule(
                    component,
                    symbols,
                    component_bonds,
                    enabled=generate_identifiers,
                )
                if identity_cache is not None:
                    identity_cache[identity_key] = identity
        else:
            identity = _DISABLED_IDENTITY
        molecule = MoleculeRecord(
            atom_ids=component,
            formula=formula,
            species_key=species_key,
            bond_signature=bond_signature,
            atom_count=len(component),
            bond_count=len(component_bonds),
            smiles=identity.smiles,
            inchi=identity.inchi,
            inchikey=identity.inchikey,
            status=identity.status,
            message=identity.message,
            instance_id=f"T{frame.timestep}:M{len(molecules):08d}",
            structure_id=topology_registry.structure_id(species_key, formula),
            atom_elements=atom_elements,
            connections=tuple(sorted(component_bonds)),
        )
        molecule_index = len(molecules)
        molecules.append(molecule)
        for atom_id in component:
            atom_to_molecule[atom_id] = molecule_index

    return FrameState(
        timestep=frame.timestep,
        bonds=bonds,
        transition_bonds=transition_bonds,
        molecules=tuple(molecules),
        atom_to_molecule=atom_to_molecule,
        symbols=symbols,
    )
