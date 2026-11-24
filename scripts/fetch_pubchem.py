"""Fetch PubChem records for the ligand list and write the local cache.

The cache is small and is committed so the screen can be rerun without
network access. Run this script again to refresh it.
"""

from __future__ import annotations

import argparse
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from druglikeness.pubchem import (  # noqa: E402
    cache_document,
    fetch_compounds,
    load_library,
    write_cache,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--library",
        type=Path,
        default=ROOT / "data" / "ligands.csv",
        help="CSV with compound_name, pubchem_cid and compound_class.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data" / "pubchem_compounds.json",
        help="JSON cache to write.",
    )
    parser.add_argument(
        "--retrieved-on",
        default=date.today().isoformat(),
        help="Retrieval date recorded in the cache (YYYY-MM-DD).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    queries = load_library(args.library)
    compounds, url = fetch_compounds(queries)
    document = cache_document(compounds, request_url=url, retrieved_on=args.retrieved_on)
    write_cache(args.output, document)
    print(f"Wrote {len(compounds)} compounds to {args.output}")
    print(f"Source: {url}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
