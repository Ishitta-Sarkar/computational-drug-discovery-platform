import json
from pathlib import Path

import pytest

from druglikeness.pubchem import load_queries, parse_property_table
from druglikeness.report import format_summary
from druglikeness.screen import screen_records


def _screen_archived(repo_root: Path):
    payload = json.loads(
        (repo_root / "results" / "pubchem_retrieval.json").read_text(encoding="utf-8")
    )
    queries = load_queries(repo_root / "data" / "ligands.csv")
    records = parse_property_table(payload["raw"], queries)
    return payload, screen_records(records)


def test_archived_structures_reproduce_the_table(repo_root: Path):
    _payload, results = _screen_archived(repo_root)
    by_name = {result.name: result for result in results}

    rutin = by_name["Rutin"]
    assert rutin.lipinski_violations == ("molecular_weight", "hbd", "hba")
    assert "logp" not in rutin.lipinski_violations
    assert rutin.veber_violations == ("tpsa",)
    assert rutin.descriptors.hbd == 10
    assert rutin.descriptors.hba == 16
    assert rutin.descriptors.molecular_weight == pytest.approx(610.521, abs=1e-3)
    assert rutin.descriptors.crippen_logp == pytest.approx(-1.6871, abs=1e-4)
    assert rutin.descriptors.rotatable_bonds == 6
    assert not rutin.retained

    quercetin = by_name["Quercetin"]
    assert quercetin.lipinski_pass
    assert quercetin.veber_pass
    assert not quercetin.qed_pass
    assert quercetin.descriptors.qed == pytest.approx(0.4341904403, abs=1e-6)

    egcg = by_name["Epigallocatechin gallate"]
    assert egcg.lipinski_violations == ("hbd", "hba")
    assert egcg.veber_violations == ("tpsa",)
    assert egcg.descriptors.rotatable_bonds == 4
    assert egcg.descriptors.hbd == 8
    assert egcg.descriptors.hba == 11

    retained = [result.name for result in results if result.retained]
    assert retained == ["Curcumin", "Resveratrol", "Luteolin"]
    assert all(0.0 <= result.descriptors.qed <= 1.0 for result in results)


def test_summary_states_the_rutin_lipinski_tension(repo_root: Path):
    payload, results = _screen_archived(repo_root)
    text = format_summary(
        results,
        retrieved_utc=payload["retrieved_utc"],
        request_url=payload["request_url"],
    )
    assert "Rutin fails Lipinski" in text
    assert "610.52" in text
    assert "-1.69" in text
    assert "donors are 10" in text
    assert "acceptors are 16" in text
    assert "Quercetin passes Lipinski" in text
    assert "disaccharide" in text
    assert "docking_score" in text
    assert "Epigallocatechin gallate" in text
    assert "0.434" in text
    assert "0.511" in text
    on_disk = (repo_root / "results" / "summary.txt").read_text(encoding="utf-8").rstrip("\n")
    assert on_disk == text


def test_readme_reports_the_computed_rutin_and_quercetin_numbers(repo_root: Path):
    readme = (repo_root / "README.md").read_text(encoding="utf-8")
    assert "610.52" in readme
    assert "0.434" in readme
    assert "0.140" in readme
    assert "results/druglikeness_screen.png" in readme
    assert "https://orcid.org/0009-0008-1028-6010" in readme
    prose = readme.casefold().replace("computational-drug-discovery-platform", "")
    assert "platform" not in prose
