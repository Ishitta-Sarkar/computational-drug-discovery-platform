"""The figure is written from screen rows and is a PNG."""

from pathlib import Path

from druglikeness.plot import plot_screen
from druglikeness.pubchem import load_cache
from druglikeness.screen import screen_compounds

ROOT = Path(__file__).resolve().parents[1]


def test_plot_writes_a_png(tmp_path: Path) -> None:
    _document, compounds = load_cache(ROOT / "data" / "pubchem_compounds.json")
    rows = screen_compounds(compounds)
    path = tmp_path / "figures" / "druglikeness_overview.png"
    plot_screen(rows, path)
    payload = path.read_bytes()
    assert payload.startswith(b"\x89PNG\r\n\x1a\n")
    assert len(payload) > 10_000
