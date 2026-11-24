# Changelog

## [0.2.0] - 2026-11-24

### Changed

- Descriptors used for filtering are calculated with RDKit from SMILES retrieved from PubChem by CID. The hand-typed molecular-weight, logP, donor and acceptor columns were removed.
- The screen is Lipinski (zero violations), Veber, and QED with a stated cutoff of 0.5.
- Docking scores were removed. The old `docking_score` column was one number per ligand, including −9.0 for rutin, −8.1 for quercetin and −7.8 for curcumin, with no target and no docking run. The quercetin and curcumin figures are the same ones listed as main-protease affinities in `covid-mpro-phytochemical-docking`. They are not results of this repository.
- The example target table (ACE2, PPARG, APP, SNCA) was removed. Those names were not structures, poses or assays.
- The write-up describes a ligand drug-likeness screening workflow.

### Added

- Pinned `requirements.txt`, pytest, and GitHub Actions (ruff and pytest).
- `CITATION.cff` and this changelog.
- Results table, text summary and figure written by `main.py`.

## [0.1.0] - 2026-06-16

### Added

- Initial ligand and target tables, a Lipinski filter on hand-entered descriptors, and ranking by an unattached docking score.
