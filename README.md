# Computational Drug Discovery Platform

A Python-based framework developed to explore fundamental computational drug discovery workflows, including target information analysis, ligand screening, docking-score ranking, and rule-based candidate filtering.

## Project Overview

Computational drug discovery plays an important role in modern pharmaceutical and biomedical research by helping researchers screen, prioritize, and analyze potential therapeutic compounds before experimental validation.

This project implements a simplified drug discovery workflow using structured ligand and target datasets. It focuses on candidate ranking based on docking scores and basic drug-likeness filtering using Lipinski-style molecular property criteria.

## Implemented Features

- Ligand library handling
- Target protein dataset handling
- Docking-score based ligand ranking
- Basic Lipinski rule-based filtering
- Candidate summary generation
- CSV-based computational workflow

## Repository Structure

```text
computational-drug-discovery-platform/

README.md
main.py
LICENSE
.gitignore

data/
    ligand_library.csv
    target_protein.csv

src/
    drug_discovery_tools.py
```

## Example Workflow

```text
Target Protein Data
        ↓
Ligand Library
        ↓
Docking Score Ranking
        ↓
Drug-Likeness Filtering
        ↓
Candidate Selection
```

## Example Compounds Included

- Quercetin
- Curcumin
- Resveratrol
- Epigallocatechin Gallate
- Luteolin
- Rutin

## Example Target Areas

- Diabetes
- Alzheimer Disease
- Parkinson Disease
- Cardiovascular Disease

## Technologies

- Python 3
- CSV Data Processing
- Object-Oriented Programming
- Computational Biology
- Computational Drug Discovery

## Scientific Context

Molecular docking and virtual screening are widely used computational approaches for prioritizing compounds based on predicted binding affinity and molecular properties. This project demonstrates a simplified version of such a workflow for learning and portfolio development.

The results from this type of computational workflow should be interpreted as preliminary and require experimental validation.

## Future Development

Planned extensions include:

- Integration with real docking output files
- Protein-ligand interaction summaries
- RDKit-based molecular descriptor calculation
- ADMET property filtering
- Visualization of ligand ranking results
- Machine learning-based bioactivity prediction
- Molecular docking workflow documentation

## Author

**Ishitta Sarkar**

B.Tech Biotechnology

Areas of Interest:

- Bioinformatics
- Computational Biology
- Molecular Docking
- Drug Discovery
- Precision Medicine
- Biomedical Data Science
