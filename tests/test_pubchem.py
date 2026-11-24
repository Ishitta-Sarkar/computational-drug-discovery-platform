from pathlib import Path

import pytest

from druglikeness.pubchem import (
    LigandQuery,
    build_url,
    load_queries,
    parse_property_table,
)


def _raw_row(**overrides):
    row = {
        "CID": 5280805,
        "Title": "Rutin",
        "IUPACName": "example",
        "SMILES": "CCO",
        "MolecularFormula": "C2H6O",
        "MolecularWeight": "46.07",
        "XLogP": -0.1,
        "TPSA": 20.2,
        "HBondDonorCount": 1,
        "HBondAcceptorCount": 1,
        "RotatableBondCount": 0,
    }
    row.update(overrides)
    return {"PropertyTable": {"Properties": [row]}}


def test_ligand_table_has_cids_and_no_docking_score(repo_root: Path):
    queries = load_queries(repo_root / "data" / "ligands.csv")
    assert [query.pubchem_cid for query in queries] == [
        5280343,
        969516,
        445154,
        65064,
        5280445,
        5280805,
    ]
    header = (repo_root / "data" / "ligands.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "docking_score" not in header
    assert not (repo_root / "data" / "target_protein.csv").exists()
    assert not (repo_root / "data" / "ligand_library.csv").exists()


def test_build_url_requests_smiles_and_reference_properties():
    url = build_url([5280805, 5280343])
    assert url.startswith("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/5280805,5280343/")
    assert "SMILES" in url
    assert "XLogP" in url


def test_parser_reads_pubchem_types_and_checks_the_title():
    records = parse_property_table(_raw_row(), [LigandQuery("Rutin", 5280805)])
    assert records[0].smiles == "CCO"
    assert records[0].reference.molecular_weight == pytest.approx(46.07)
    assert records[0].reference.xlogp == pytest.approx(-0.1)
    assert records[0].reference.hbd == 1


def test_parser_rejects_a_title_that_does_not_match_the_query():
    with pytest.raises(ValueError, match="titled"):
        parse_property_table(_raw_row(Title="Quercetin"), [LigandQuery("Rutin", 5280805)])


def test_parser_rejects_a_missing_cid():
    with pytest.raises(ValueError, match="no CID"):
        parse_property_table(_raw_row(), [LigandQuery("Quercetin", 5280343)])


def test_load_queries_rejects_a_docking_column(tmp_path: Path):
    path = tmp_path / "ligands.csv"
    path.write_text(
        "compound_name,pubchem_cid,docking_score\nRutin,5280805,-9.0\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="docking_score"):
        load_queries(path)
