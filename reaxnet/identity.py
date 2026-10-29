from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from reaxnet.models import Bond, BondKey


@dataclass(frozen=True)
class IdentityResult:
    smiles: str = ""
    inchi: str = ""
    inchikey: str = ""
    status: str = "unavailable"
    message: str = ""


def identify_molecule(
    atom_ids: tuple[int, ...],
    symbols: Mapping[int, str],
    bonds: Mapping[BondKey, Bond],
    *,
    enabled: bool = True,
) -> IdentityResult:
    if not enabled:
        return IdentityResult(status="disabled", message="Identifier generation was disabled")

    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem import inchi as rd_inchi
    except ImportError:
        return IdentityResult(
            status="unavailable",
            message="RDKit is not installed; graph analysis completed without SMILES/InChI",
        )

    def build_molecule(use_bond_orders: bool):
        editable = Chem.RWMol()
        atom_indexes: dict[int, int] = {}
        for atom_id in atom_ids:
            atom = Chem.Atom(symbols[atom_id])
            atom.SetNoImplicit(True)
            atom_indexes[atom_id] = editable.AddAtom(atom)

        for key, bond in sorted(bonds.items()):
            if key[0] not in atom_indexes or key[1] not in atom_indexes:
                continue
            bond_type = Chem.BondType.SINGLE
            if use_bond_orders:
                if bond.effective_order >= 2.5:
                    bond_type = Chem.BondType.TRIPLE
                elif bond.effective_order >= 1.5:
                    bond_type = Chem.BondType.DOUBLE
            editable.AddBond(atom_indexes[key[0]], atom_indexes[key[1]], bond_type)

        molecule = editable.GetMol()
        Chem.SanitizeMol(molecule)
        return molecule

    fallback_message = ""
    with rdBase.BlockLogs():
        try:
            molecule = build_molecule(use_bond_orders=True)
        except Exception as primary_error:
            return IdentityResult(status="failed", message="RDKit could not build the molecule: " +
                                  " ".join(str(primary_error).split()))

        try:
            display_molecule = Chem.RemoveHs(molecule)
            smiles = Chem.MolToSmiles(display_molecule, canonical=True)
        except Exception as exc:
            return IdentityResult(
                status="failed",
                message=f"RDKit could not generate SMILES: {' '.join(str(exc).split())}",
            )

        try:
            inchi = rd_inchi.MolToInchi(molecule)
            inchikey = rd_inchi.InchiToInchiKey(inchi) if inchi else ""
        except Exception as exc:
            inchi = ""
            inchikey = ""
            message = " ".join(str(exc).split())
            return IdentityResult(
                smiles=smiles,
                status="partial",
                message=f"SMILES generated, but InChI generation failed: {message}",
            )

    if not inchi or not inchikey:
        return IdentityResult(
            smiles=smiles,
            inchi=inchi,
            inchikey=inchikey,
            status="partial",
            message=fallback_message or "SMILES generated, but InChI or InChIKey is empty",
        )
    return IdentityResult(
        smiles=smiles,
        inchi=inchi,
        inchikey=inchikey,
        status="fallback" if fallback_message else "ok",
        message=fallback_message,
    )

@dataclass(frozen=True)
class CandidateIdentifiers:
    canonical_smiles: str = ''
    inchi: str = ''
    inchikey: str = ''
    smiles_status: str = 'unavailable'
    inchi_status: str = 'unavailable'
    validation_status: str = 'unchecked'
    message: str = ''


def identify_structure_candidate(symbols, edges, bond_orders, charges, radicals):
    """Validate an explicit candidate without changing its atoms or electron state.

    Indices in edges refer to symbols. No BO rounding or single-bond fallback.
    """
    from collections import Counter
    try:
        from rdkit import Chem, rdBase
    except ImportError:
        return CandidateIdentifiers(message='RDKit unavailable; candidate graph retained')
    with rdBase.BlockLogs():
        try:
            if not (len(symbols) == len(charges) == len(radicals)) or len(edges) != len(bond_orders):
                raise ValueError('candidate array lengths disagree')
            editable = Chem.RWMol()
            for symbol, charge, radical in zip(symbols, charges, radicals):
                atom = Chem.Atom(symbol)
                atom.SetNoImplicit(True)
                atom.SetFormalCharge(charge)
                atom.SetNumRadicalElectrons(radical)
                editable.AddAtom(atom)
            types = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
            for (a, b), order in zip(edges, bond_orders):
                editable.AddBond(a, b, types[order])
            molecule = editable.GetMol()
            Chem.SanitizeMol(molecule)
            if tuple(a.GetSymbol() for a in molecule.GetAtoms()) != tuple(symbols):
                raise ValueError('element sequence changed')
            if tuple(a.GetFormalCharge() for a in molecule.GetAtoms()) != tuple(charges):
                raise ValueError('formal charges changed during sanitization')
            if tuple(a.GetNumRadicalElectrons() for a in molecule.GetAtoms()) != tuple(radicals):
                raise ValueError('radicals changed during sanitization')
            actual_edges = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in molecule.GetBonds()}
            if actual_edges != set(edges):
                raise ValueError('adjacency changed during sanitization')
            if any(a.GetNumImplicitHs() or a.GetNumExplicitHs() for a in molecule.GetAtoms()):
                raise ValueError('sanitization introduced unobserved hydrogen counts')
        except Exception as exc:
            return CandidateIdentifiers(validation_status='failed', smiles_status='not_attempted',
                                        inchi_status='not_attempted', message=str(exc))
        try:
            display = Chem.RemoveHs(molecule)
            smiles = Chem.MolToSmiles(display, canonical=True, isomericSmiles=True)
            restored = Chem.MolFromSmiles(smiles)
            if restored is None:
                raise ValueError('SMILES cannot be read back')
            restored = Chem.AddHs(restored)
            if Counter(a.GetSymbol() for a in restored.GetAtoms()) != Counter(symbols):
                raise ValueError('SMILES round trip changed element or hydrogen counts')
            original_key = Chem.MolToSmiles(molecule, canonical=True, allHsExplicit=True)
            restored_key = Chem.MolToSmiles(restored, canonical=True, allHsExplicit=True)
            if original_key != restored_key:
                raise ValueError('SMILES round trip changed the explicit molecular graph')
        except Exception as exc:
            return CandidateIdentifiers(validation_status='passed', smiles_status='failed',
                                        inchi_status='not_attempted', message=str(exc))
        try:
            from rdkit.Chem import inchi as rd_inchi
            inchi = rd_inchi.MolToInchi(molecule)
            inchikey = rd_inchi.InchiToInchiKey(inchi) if inchi else ''
            if not inchi or not inchikey:
                raise ValueError('InChI or InChIKey is empty')
        except Exception as exc:
            return CandidateIdentifiers(canonical_smiles=smiles, smiles_status='ok', inchi_status='failed',
                                        validation_status='passed', message=str(exc))
    return CandidateIdentifiers(smiles, inchi, inchikey, 'ok', 'ok', 'passed')
