"""RDKit descriptors used by the drug-likeness rules.

The numbers are calculated from a SMILES string. They are not typed in.
Lipinski hydrogen-bond counts use RDKit's Lipinski definitions. logP is the
Wildman–Crippen logP (RDKit MolLogP), which is not an experimental logP and
is not the Moriguchi logP cited in some statements of the rule of five.
QED uses RDKit's default mean descriptor weights.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import QED, Descriptors, inchi, rdMolDescriptors


@dataclass(frozen=True)
class MolecularDescriptors:
    """Descriptors for one parsed molecule."""

    molecular_weight: float
    crippen_logp: float
    hbd_lipinski: int
    hba_lipinski: int
    tpsa: float
    rotatable_bonds: int
    qed: float
    inchikey: str


def mol_from_smiles(smiles: str) -> Chem.Mol:
    """Parse SMILES. Raise ValueError when RDKit cannot build a molecule."""
    if not smiles or not smiles.strip():
        raise ValueError("SMILES string is empty.")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    return mol


def descriptors_from_mol(mol: Chem.Mol) -> MolecularDescriptors:
    """Compute the descriptor set used by the filters."""
    qed_value = float(QED.qed(mol))
    if not 0.0 <= qed_value <= 1.0:
        raise ValueError(f"RDKit returned a QED outside [0, 1]: {qed_value}")
    return MolecularDescriptors(
        molecular_weight=float(Descriptors.MolWt(mol)),
        crippen_logp=float(Descriptors.MolLogP(mol)),
        hbd_lipinski=int(rdMolDescriptors.CalcNumLipinskiHBD(mol)),
        hba_lipinski=int(rdMolDescriptors.CalcNumLipinskiHBA(mol)),
        tpsa=float(Descriptors.TPSA(mol)),
        rotatable_bonds=int(Descriptors.NumRotatableBonds(mol)),
        qed=qed_value,
        inchikey=inchi.MolToInchiKey(mol),
    )


def descriptors_from_smiles(smiles: str) -> MolecularDescriptors:
    """Parse SMILES and compute descriptors."""
    return descriptors_from_mol(mol_from_smiles(smiles))
