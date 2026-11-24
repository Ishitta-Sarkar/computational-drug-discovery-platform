import json
from pathlib import Path

from druglikeness.descriptors import compute_descriptors
from druglikeness.pubchem import (
    CompoundRecord,
    LigandQuery,
    ReferenceProperties,
    Retrieval,
    fetch_records,
    load_queries,
)
from druglikeness.run import run
from druglikeness.screen import screen_records


def _ethanol_retrieval(queries: list[LigandQuery]) -> Retrieval:
    assert len(queries) == 1
    query = queries[0]
    reference = ReferenceProperties(46.07, -0.1, 20.2, 1, 1, 0)
    record = CompoundRecord(
        name=query.name,
        pubchem_cid=query.pubchem_cid,
        pubchem_title=query.name,
        iupac_name="ethanol",
        smiles="CCO",
        molecular_formula="C2H6O",
        reference=reference,
    )
    return Retrieval(
        retrieved_utc="2026-11-24T00:00:00+00:00",
        request_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/example/property/SMILES/JSON",
        raw={"PropertyTable": {"Properties": []}},
        records=(record,),
    )


def test_run_writes_table_summary_and_figure(tmp_path: Path):
    data = tmp_path / "data"
    data.mkdir()
    (data / "ligands.csv").write_text("compound_name,pubchem_cid\nEthanol,702\n", encoding="utf-8")
    run(tmp_path, fetcher=_ethanol_retrieval)
    table = (tmp_path / "results" / "druglikeness_table.csv").read_text(encoding="utf-8")
    assert "docking_score" not in table
    assert "Ethanol" in table
    summary = (tmp_path / "results" / "summary.txt").read_text(encoding="utf-8")
    assert "Ethanol" in summary
    assert "does not dock" in summary
    figure = tmp_path / "results" / "druglikeness_screen.png"
    assert figure.read_bytes().startswith(b"\x89PNG")
    stored = json.loads(
        (tmp_path / "results" / "pubchem_retrieval.json").read_text(encoding="utf-8")
    )
    assert stored["records"][0]["smiles"] == "CCO"


def test_live_pubchem_fetch_matches_the_ligand_table(repo_root: Path):
    queries = load_queries(repo_root / "data" / "ligands.csv")
    retrieval = fetch_records(queries)
    assert [record.pubchem_cid for record in retrieval.records] == [
        query.pubchem_cid for query in queries
    ]
    assert "pubchem.ncbi.nlm.nih.gov" in retrieval.request_url
    results = screen_records(retrieval.records)
    rutin = next(result for result in results if result.name == "Rutin")
    assert rutin.molecular_formula == "C27H30O16"
    assert rutin.lipinski_violations == ("molecular_weight", "hbd", "hba")
    quercetin = next(result for result in results if result.name == "Quercetin")
    assert compute_descriptors(quercetin.smiles, quercetin.molecular_formula).hbd == 5
