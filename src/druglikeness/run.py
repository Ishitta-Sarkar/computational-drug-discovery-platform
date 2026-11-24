"""Run the screen and write the table, summary and figure."""

import json
from pathlib import Path

from druglikeness.plot import plot_screen
from druglikeness.pubchem import Retrieval, fetch_records, load_queries
from druglikeness.report import format_summary
from druglikeness.screen import screen_records, write_table


def run(root: Path, fetcher=None) -> None:
    queries = load_queries(root / "data" / "ligands.csv")
    retrieval = (fetcher or fetch_records)(queries)
    results = screen_records(retrieval.records)
    output = root / "results"
    output.mkdir(parents=True, exist_ok=True)
    _write_retrieval(output / "pubchem_retrieval.json", retrieval)
    write_table(output / "druglikeness_table.csv", results)
    summary = format_summary(
        results,
        retrieved_utc=retrieval.retrieved_utc,
        request_url=retrieval.request_url,
    )
    (output / "summary.txt").write_text(summary + "\n", encoding="utf-8")
    plot_screen(results, output / "druglikeness_screen.png")
    print(summary)


def _write_retrieval(path: Path, retrieval: Retrieval) -> None:
    payload = {
        "source": "PubChem PUG REST",
        "smiles_field": "SMILES",
        "retrieved_utc": retrieval.retrieved_utc,
        "request_url": retrieval.request_url,
        "records": [
            {
                "compound_name": record.name,
                "pubchem_cid": record.pubchem_cid,
                "pubchem_title": record.pubchem_title,
                "iupac_name": record.iupac_name,
                "smiles": record.smiles,
                "molecular_formula": record.molecular_formula,
                "reference": {
                    "molecular_weight": record.reference.molecular_weight,
                    "xlogp": record.reference.xlogp,
                    "tpsa": record.reference.tpsa,
                    "hbd": record.reference.hbd,
                    "hba": record.reference.hba,
                    "rotatable_bonds": record.reference.rotatable_bonds,
                },
            }
            for record in retrieval.records
        ],
        "raw": retrieval.raw,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
