from rdkit import Chem

from druglikeness.descriptors import compute_descriptors, veber_rotatable_bond_count


def test_amide_is_excluded_and_ester_is_counted():
    n_methylacetamide = Chem.MolFromSmiles("CC(=O)NC")
    methyl_acetate = Chem.MolFromSmiles("CC(=O)OC")
    assert veber_rotatable_bond_count(n_methylacetamide) == 0
    assert veber_rotatable_bond_count(methyl_acetate) == 1


def test_explicit_hydrogens_do_not_change_the_count():
    mol = Chem.MolFromSmiles("CCc1ccccc1")
    assert veber_rotatable_bond_count(mol) == veber_rotatable_bond_count(Chem.AddHs(mol))


def test_formula_mismatch_is_rejected():
    try:
        compute_descriptors("CCO", expected_formula="C3H8")
    except ValueError as exc:
        assert "Formula mismatch" in str(exc)
    else:
        raise AssertionError("expected a formula mismatch")


def test_unparseable_smiles_is_rejected():
    try:
        compute_descriptors("not-a-molecule")
    except ValueError as exc:
        assert "could not parse" in str(exc)
    else:
        raise AssertionError("expected a parse error")
