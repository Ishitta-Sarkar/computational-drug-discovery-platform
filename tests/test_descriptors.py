"""Descriptor calculations on structures that are not the project ligand set."""

import pytest
from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from druglikeness.descriptors import descriptors_from_smiles, mol_from_smiles
from druglikeness.filters import lipinski_assessment, pains_descriptions, veber_assessment

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
ETHANE = "CC"
CATECHOL = "c1ccc(O)c(O)c1"


def test_empty_and_invalid_smiles_are_rejected() -> None:
    with pytest.raises(ValueError, match="empty"):
        mol_from_smiles("  ")
    with pytest.raises(ValueError, match="could not parse"):
        descriptors_from_smiles("not-a-molecule")


def test_aspirin_descriptors_match_rdkit_definitions() -> None:
    desc = descriptors_from_smiles(ASPIRIN)
    assert desc.molecular_weight == pytest.approx(180.16, abs=0.02)
    assert desc.hbd_lipinski == 1
    assert desc.hba_lipinski == 4
    assert desc.rotatable_bonds == 2
    assert 0.0 <= desc.qed <= 1.0
    assert lipinski_assessment(desc).passed
    assert veber_assessment(desc).passed
    assert pains_descriptions(mol_from_smiles(ASPIRIN)) == []


def test_ethane_has_no_pains_alert() -> None:
    assert pains_descriptions(mol_from_smiles(ETHANE)) == []


def test_catechol_matches_the_pains_catechol_filter() -> None:
    alerts = pains_descriptions(mol_from_smiles(CATECHOL))
    assert alerts == ["catechol_A(92)"]
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    matches = catalog.GetMatches(Chem.MolFromSmiles(CATECHOL))
    descriptions = [entry.GetDescription() for entry in matches]
    assert descriptions == alerts
