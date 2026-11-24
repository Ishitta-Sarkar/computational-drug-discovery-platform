"""Command-line screen against the committed cache."""

from pathlib import Path

import main

from druglikeness.screen import RESULT_FIELDS


def test_main_writes_table_figure_and_run_info(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    results = tmp_path / "screen_results.csv"
    figure = tmp_path / "druglikeness_overview.png"
    run_info = tmp_path / "run_info.json"
    exit_code = main.main(
        [
            "--library",
            str(root / "data" / "ligands.csv"),
            "--cache",
            str(root / "data" / "pubchem_compounds.json"),
            "--results",
            str(results),
            "--figure",
            str(figure),
            "--run-info",
            str(run_info),
        ]
    )
    assert exit_code == 0
    header = results.read_text(encoding="utf-8").splitlines()[0]
    assert header.split(",") == RESULT_FIELDS
    assert "docking_score" not in results.read_text(encoding="utf-8")
    assert figure.read_bytes().startswith(b"\x89PNG")
    text = run_info.read_text(encoding="utf-8")
    assert "2026.09.2" in text
    assert "pubchem" in text
