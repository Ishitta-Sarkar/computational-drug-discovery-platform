"""Screen of the committed PubChem cache.

Expected values are RDKit results for these SMILES at the pinned RDKit
version, not the hand-typed numbers from the previous ligand table.
"""

from pathlib import Path

import pytest

from druglikeness.pubchem import load_cache
from druglikeness.reporting import format_report
from druglikeness.screen import screen_compounds, write_results

ROOT = Path(__file__).resolve().parents[1]


def _rows() -> list[dict[str, object]]:
    _document, compounds = load_cache(ROOT / "data" / "pubchem_compounds.json")
    return screen_compounds(compounds)


def _by_name(rows: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    return {str(row["compound_name"]): row for row in rows}


def test_cache_screen_has_no_docking_score() -> None:
    rows = _rows()
    assert len(rows) == 6
    for row in rows:
        assert "docking_score" not in row
        assert "kcal" not in row
        assert row["rdkit_version"]


def test_rutin_fails_lipinski_and_veber_tpsa_while_quercetin_passes() -> None:
    rows = _by_name(_rows())
    rutin = rows["Rutin"]
    quercetin = rows["Quercetin"]

    assert rutin["molecular_weight"] == pytest.approx(610.521, abs=0.01)
    assert rutin["hbd_lipinski"] == 10
    assert rutin["hba_lipinski"] == 16
    assert float(rutin["crippen_logp"]) < 0
    assert rutin["lipinski_violations"] == "molecular_weight;hbd_lipinski;hba_lipinski"
    assert rutin["lipinski_pass_strict"] is False
    assert rutin["lipinski_alert"] is True
    assert rutin["rotatable_bonds"] <= 10
    assert rutin["veber_violations"] == "tpsa"
    assert float(rutin["tpsa"]) > 140

    assert quercetin["hbd_lipinski"] == 5
    assert quercetin["hba_lipinski"] == 7
    assert quercetin["molecular_weight"] == pytest.approx(302.238, abs=0.01)
    assert quercetin["lipinski_pass_strict"] is True
    assert quercetin["lipinski_alert"] is False
    assert quercetin["veber_pass"] is True
    assert float(rutin["qed"]) < float(quercetin["qed"])


def test_egcg_fails_hydrogen_bonding_rules_and_not_molecular_weight() -> None:
    egcg = _by_name(_rows())["Epigallocatechin gallate"]
    assert float(egcg["molecular_weight"]) < 500
    assert egcg["hbd_lipinski"] == 8
    assert egcg["hba_lipinski"] == 11
    assert egcg["lipinski_violations"] == "hbd_lipinski;hba_lipinski"
    assert egcg["veber_violations"] == "tpsa"
    assert egcg["pains_alerts"] == "catechol_A(92)"


def test_pass_and_pains_sets() -> None:
    rows = _by_name(_rows())
    strict = {name for name, row in rows.items() if row["lipinski_pass_strict"]}
    assert strict == {"Quercetin", "Curcumin", "Resveratrol", "Luteolin"}
    pains = {name for name, row in rows.items() if row["pains_flag"]}
    assert pains == {"Quercetin", "Epigallocatechin gallate", "Luteolin", "Rutin"}
    assert rows["Curcumin"]["pains_alerts"] == ""
    assert rows["Resveratrol"]["pains_alerts"] == ""
    assert all(row["veber_pass"] for name, row in rows.items() if name in strict)
    assert float(rows["Resveratrol"]["qed"]) == max(float(row["qed"]) for row in rows.values())
    assert float(rows["Rutin"]["qed"]) == min(float(row["qed"]) for row in rows.values())


def test_every_logp_in_this_set_is_inside_the_lipinski_limit() -> None:
    rows = _rows()
    assert all(float(row["crippen_logp"]) <= 5 for row in rows)
    assert all("crippen_logp" not in str(row["lipinski_violations"]) for row in rows)


def test_committed_results_match_a_fresh_screen(tmp_path: Path) -> None:
    output = tmp_path / "screen_results.csv"
    write_results(output, _rows())
    committed = ROOT / "results" / "screen_results.csv"
    assert output.read_text(encoding="utf-8") == committed.read_text(encoding="utf-8")


def test_report_states_the_rutin_result_and_does_not_rank_by_a_score() -> None:
    _document, compounds = load_cache(ROOT / "data" / "pubchem_compounds.json")
    report = format_report(screen_compounds(compounds), retrieved_on="2026-11-24")
    assert "Rutin Lipinski: molecular weight; H-bond donors; H-bond acceptors." in report
    assert "Quercetin Lipinski: pass." in report
    assert "docking" in report
    assert "-9.0" in report
    table = report.split("Compound", maxsplit=1)[1].split("Lipinski pass", maxsplit=1)[0]
    assert "kcal/mol" not in table
