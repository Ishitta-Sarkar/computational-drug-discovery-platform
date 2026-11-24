"""Run the ligand drug-likeness screen from a fresh clone."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from druglikeness.run import run


def main() -> None:
    run(ROOT)


if __name__ == "__main__":
    main()
