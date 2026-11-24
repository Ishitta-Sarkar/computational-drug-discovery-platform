"""Run the drug-likeness screen and write the table and figure.

Descriptors are computed with RDKit from the SMILES in the PubChem cache.
The cache is produced by ``scripts/fetch_pubchem.py``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from druglikeness.plot import plot_screen  # noqa: E402
from druglikeness.pubchem import (  # noqa: E402
    assert_cache_matches_library,
    load_cache,
    load_library,
)
from druglikeness.reporting import format_report  # noqa: E402
from druglikeness.screen import screen_compounds, write_results  # noqa: E402


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path, default=ROOT / "data" / "ligands.csv")
    parser.add_argument("--cache", type=Path, default=ROOT / "data" / "pubchem_compounds.json")
    parser.add_argument("--results", type=Path, default=ROOT / "results" / "screen_results.csv")
    parser.add_argument(
        "--figure",
        type=Path,
        default=ROOT / "figures" / "druglikeness_overview.png",
    )
    parser.add_argument("--run-info", type=Path, default=ROOT / "results" / "run_info.json")
    return parser.parse_args(argv)


def _relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    queries = load_library(args.library)
    document, compounds = load_cache(args.cache)
    assert_cache_matches_library(compounds, queries)
    rows = screen_compounds(compounds)
    write_results(args.results, rows)
    plot_screen(rows, args.figure)
    run_info = {
        "rdkit_version": rows[0]["rdkit_version"],
        "n_compounds": len(rows),
        "cache_retrieved_on": document["retrieved_on"],
        "cache_request_url": document["request_url"],
        "library": _relative(args.library),
        "cache": _relative(args.cache),
        "results_csv": _relative(args.results),
        "figure": _relative(args.figure),
    }
    args.run_info.parent.mkdir(parents=True, exist_ok=True)
    args.run_info.write_text(json.dumps(run_info, indent=2) + "\n", encoding="utf-8")
    print(format_report(rows, retrieved_on=str(document["retrieved_on"])))
    print()
    print(f"Wrote {args.results}")
    print(f"Wrote {args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
