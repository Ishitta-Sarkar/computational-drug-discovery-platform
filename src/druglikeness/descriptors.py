"""RDKit descriptors used by the drug-likeness filters."""

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.QED import qed


@dataclass(frozen=True)
class DescriptorSet:
    """Descriptors calculated for one structure.

    ``crippen_logp`` is RDKit's Wildman–Crippen logP (``Crippen.MolLogP``),
    the same lipophilicity estimate QED uses. It is not PubChem XLogP.
    ``rotatable_bonds`` follows Veber et al. (see ``veber_rotatable_bond_count``)
    and is not the strict count RDKit uses inside QED.
    """

    molecular_weight: float
    crippen_logp: float
    hbd: int
    hba: int
    tpsa: float
    rotatable_bonds: int
    qed: float


def molecular_formula(mol: Chem.Mol) -> str:
    return rdMolDescriptors.CalcMolFormula(mol)


def veber_rotatable_bond_count(mol: Chem.Mol) -> int:
    """Count rotatable bonds as defined by Veber et al., J Med Chem 2002.

    A bond counts when it is a single bond, is not in a ring, and both atoms
    are non-terminal heavy atoms. An atom is terminal when it has fewer than
    two heavy-atom neighbors. Amide C–N bonds are excluded because of their
    rotational barrier. Ester bonds are counted.

    RDKit's ``Lipinski.NumRotatableBonds`` uses the strict count, which also
    drops esters. That count is what QED uses internally. It is not used here.
    """

    count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
            continue
        begin = bond.GetBeginAtom()
        end = bond.GetEndAtom()
        if begin.GetAtomicNum() == 1 or end.GetAtomicNum() == 1:
            continue
        if _heavy_degree(begin) < 2 or _heavy_degree(end) < 2:
            continue
        if _is_amide_cn(bond):
            continue
        count += 1
    return count


def compute_descriptors(smiles: str, expected_formula: str | None = None) -> DescriptorSet:
    """Calculate the descriptor set for one SMILES string.

    When ``expected_formula`` is given, the RDKit formula must match it.
    That check catches a CID that was paired with the wrong structure.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    formula = molecular_formula(mol)
    if expected_formula is not None and formula != expected_formula:
        raise ValueError(
            f"Formula mismatch for SMILES {smiles}: RDKit has {formula}, "
            f"expected {expected_formula}."
        )
    return DescriptorSet(
        molecular_weight=float(Descriptors.MolWt(mol)),
        crippen_logp=float(Crippen.MolLogP(mol)),
        hbd=int(Lipinski.NumHDonors(mol)),
        hba=int(Lipinski.NumHAcceptors(mol)),
        tpsa=float(Descriptors.TPSA(mol)),
        rotatable_bonds=veber_rotatable_bond_count(mol),
        qed=float(qed(mol)),
    )


def _heavy_degree(atom: Chem.Atom) -> int:
    return sum(1 for neighbour in atom.GetNeighbors() if neighbour.GetAtomicNum() > 1)


def _is_amide_cn(bond: Chem.Bond) -> bool:
    begin = bond.GetBeginAtom()
    end = bond.GetEndAtom()
    if {begin.GetAtomicNum(), end.GetAtomicNum()} != {6, 7}:
        return False
    carbon = begin if begin.GetAtomicNum() == 6 else end
    mol = bond.GetOwningMol()
    for neighbour in carbon.GetNeighbors():
        if neighbour.GetAtomicNum() != 8:
            continue
        carbonyl = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbour.GetIdx())
        if carbonyl is not None and carbonyl.GetBondType() == Chem.BondType.DOUBLE:
            return True
    return False
