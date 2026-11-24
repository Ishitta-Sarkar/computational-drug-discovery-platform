"""Ligand-list parsing and PubChem property-table alignment. No network."""

import json
from pathlib import Path

import pytest

from druglikeness.pubchem import (
    LigandQuery,
    assert_cache_matches_library,
    cache_document,
    load_cache,
    load_library,
    pug_property_url,
    records_from_property_table,
    write_cache,
)

ROOT = Path(__file__).resolve().parents[1]


def _query() -> LigandQuery:
    return LigandQuery("Quercetin", 5280343, "flavonol")


def _property_row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "CID": 5280343,
        "Title": "Quercetin",
        "IUPACName": "2-(3,4-dihydroxyphenyl)-3,5,7-trihydroxychromen-4-one",
        "MolecularFormula": "C15H10O7",
        "SMILES": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        "ConnectivitySMILES": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        "InChI": (
            "InChI=1S/C15H10O7/c16-7-4-10(19)12-11(5-7)22-15(14(21)13(12)20)"
            "6-1-2-8(17)9(18)3-6/h1-5,16-19,21H"
        ),
        "InChIKey": "REFJWTPEDVJJIY-UHFFFAOYSA-N",
    }
    row.update(overrides)
    return row


def test_library_loads_six_cids_and_rejects_a_docking_column(tmp_path: Path) -> None:
    queries = load_library(ROOT / "data" / "ligands.csv")
    assert [query.compound_name for query in queries] == [
        "Quercetin",
        "Curcumin",
        "Resveratrol",
        "Epigallocatechin gallate",
        "Luteolin",
        "Rutin",
    ]
    assert queries[0].pubchem_cid == 5280343
    bad = tmp_path / "ligands.csv"
    bad.write_text(
        "compound_name,pubchem_cid,compound_class,docking_score\nQuercetin,5280343,flavonol,-8.1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="docking_score"):
        load_library(bad)


def test_duplicate_cid_is_rejected(tmp_path: Path) -> None:
    path = tmp_path / "ligands.csv"
    path.write_text(
        "compound_name,pubchem_cid,compound_class\n"
        "Quercetin,5280343,flavonol\n"
        "Quercetin again,5280343,flavonol\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="more than once"):
        load_library(path)


def test_property_url_lists_every_cid() -> None:
    url = pug_property_url([5280343, 5280805])
    assert url.startswith("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/")
    assert "5280343,5280805" in url
    assert "IsomericSMILES" not in url
    assert "SMILES" in url


def test_title_mismatch_and_missing_smiles_fail() -> None:
    with pytest.raises(ValueError, match="does not match"):
        records_from_property_table([_query()], [_property_row(Title="Kaempferol")])
    with pytest.raises(ValueError, match="missing SMILES"):
        records_from_property_table([_query()], [_property_row(SMILES="")])


def test_case_only_title_difference_is_accepted() -> None:
    compounds = records_from_property_table(
        [LigandQuery("Epigallocatechin gallate", 65064, "flavan-3-ol gallate")],
        [
            {
                "CID": 65064,
                "Title": "Epigallocatechin Gallate",
                "IUPACName": "example",
                "MolecularFormula": "C22H18O11",
                "SMILES": "C",
                "ConnectivitySMILES": "C",
                "InChI": "InChI=1S/CH4/h1H4",
                "InChIKey": "VNWKTOKETHGBQD-UHFFFAOYSA-N",
            }
        ],
    )
    assert compounds[0].compound_name == "Epigallocatechin gallate"
    assert compounds[0].title == "Epigallocatechin Gallate"


def test_cache_roundtrip_and_library_check(tmp_path: Path) -> None:
    compounds = records_from_property_table([_query()], [_property_row()])
    url = pug_property_url([5280343])
    document = cache_document(compounds, request_url=url, retrieved_on="2026-11-24")
    path = tmp_path / "cache.json"
    write_cache(path, document)
    loaded_document, loaded = load_cache(path)
    assert loaded_document["retrieved_on"] == "2026-11-24"
    assert loaded[0] == compounds[0]
    assert_cache_matches_library(loaded, [_query()])
    with pytest.raises(ValueError, match="does not match"):
        assert_cache_matches_library(loaded, [LigandQuery("Rutin", 5280805, "flavonol glycoside")])
    broken = json.loads(path.read_text(encoding="utf-8"))
    broken["schema_version"] = 99
    path.write_text(json.dumps(broken), encoding="utf-8")
    with pytest.raises(ValueError, match="schema_version"):
        load_cache(path)
