# Changelog

## [0.2.0] - 2026-11-24

### Changed

- Molecular descriptors are computed with RDKit from PubChem isomeric SMILES. The ligand file stores names and PubChem CIDs. A script fetches the structures.
- The screen reports Lipinski violations (strict zero-violation pass, and the original two-or-more alert), Veber violations, QED and PAINS alerts.
- The README describes this repository as a learning project and a pre-docking filter.

### Removed

- Per-ligand docking scores. The previous table had one score per ligand and no target, so the score was not a ligand–target result. No docking run is stored in this repository.
- The example target table (ACE2, PPARG, APP, SNCA). Those proteins were not docked here.

### Added

- PubChem cache with request URL and retrieval date, checked against RDKit InChIKeys.
- Tests, Ruff lint, GitHub Actions CI, pinned requirements, CITATION.cff, and a figure written by the screen.

### Fixed

- Pytest can import `main.py` when it is launched from another working directory, which is how the CI job runs it.

## [0.1.0] - 2026-06-16

Initial ligand table with hand-typed descriptors and a single docking score per ligand.
