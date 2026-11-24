"""Run descriptors and filters for each cached PubChem compound."""

from __future__ import annotations

import csv
from pathlib import Path

import rdkit

from druglikeness.descriptors import descriptors_from_mol, mol_from_smiles
from druglikeness.filters import (
    lipinski_alert,
    lipinski_assessment,
    pains_descriptions,
    veber_assessment,
)
from druglikeness.pubchem import PubChemCompound

RESULT_FIELDS = [
    "compound_name",
    "compound_class",
    "pubchem_cid",
    "pubchem_title",
    "iupac_name",
    "molecular_formula",
    "isomeric_smiles",
    "connectivity_smiles",
    "inchi",
    "inchikey",
    "molecular_weight",
    "crippen_logp",
    "hbd_lipinski",
    "hba_lipinski",
    "tpsa",
    "rotatable_bonds",
    "qed",
    "lipinski_violations",
    "lipinski_violation_count",
    "lipinski_pass_strict",
    "lipinski_alert",
    "veber_violations",
    "veber_pass",
    "pains_alerts",
    "pains_flag",
    "rdkit_version",
]


def _bool_text(value: bool) -> str:
    return "true" if value else "false"


def screen_compound(compound: PubChemCompound) -> dict[str, object]:
    """Compute descriptors and filters for one compound.

    The InChIKey RDKit derives from the PubChem isomeric SMILES must match
    the InChIKey PubChem reported. A mismatch means the SMILES and the
    identifier are not the same structure, and the screen stops.
    """
    mol = mol_from_smiles(compound.isomeric_smiles)
    desc = descriptors_from_mol(mol)
    if desc.inchikey != compound.inchikey:
        raise ValueError(
            f"InChIKey mismatch for {compound.compound_name} "
            f"(CID {compound.pubchem_cid}): RDKit computed {desc.inchikey} "
            f"from the PubChem SMILES, but PubChem reports {compound.inchikey}."
        )
    lipinski = lipinski_assessment(desc)
    veber = veber_assessment(desc)
    pains = pains_descriptions(mol)
    return {
        "compound_name": compound.compound_name,
        "compound_class": compound.compound_class,
        "pubchem_cid": compound.pubchem_cid,
        "pubchem_title": compound.title,
        "iupac_name": compound.iupac_name,
        "molecular_formula": compound.molecular_formula,
        "isomeric_smiles": compound.isomeric_smiles,
        "connectivity_smiles": compound.connectivity_smiles,
        "inchi": compound.inchi,
        "inchikey": compound.inchikey,
        "molecular_weight": desc.molecular_weight,
        "crippen_logp": desc.crippen_logp,
        "hbd_lipinski": desc.hbd_lipinski,
        "hba_lipinski": desc.hba_lipinski,
        "tpsa": desc.tpsa,
        "rotatable_bonds": desc.rotatable_bonds,
        "qed": desc.qed,
        "lipinski_violations": ";".join(lipinski.violations),
        "lipinski_violation_count": lipinski.count,
        "lipinski_pass_strict": lipinski.passed,
        "lipinski_alert": lipinski_alert(lipinski),
        "veber_violations": ";".join(veber.violations),
        "veber_pass": veber.passed,
        "pains_alerts": ";".join(pains),
        "pains_flag": bool(pains),
        "rdkit_version": rdkit.__version__,
    }


def screen_compounds(compounds: list[PubChemCompound]) -> list[dict[str, object]]:
    """Screen compounds in the order they are given."""
    if not compounds:
        raise ValueError("No compounds to screen.")
    return [screen_compound(compound) for compound in compounds]


def _csv_row(row: dict[str, object]) -> dict[str, object]:
    formatted = dict(row)
    formatted["molecular_weight"] = f"{row['molecular_weight']:.4f}"
    formatted["crippen_logp"] = f"{row['crippen_logp']:.4f}"
    formatted["tpsa"] = f"{row['tpsa']:.4f}"
    formatted["qed"] = f"{row['qed']:.6f}"
    formatted["lipinski_pass_strict"] = _bool_text(bool(row["lipinski_pass_strict"]))
    formatted["lipinski_alert"] = _bool_text(bool(row["lipinski_alert"]))
    formatted["veber_pass"] = _bool_text(bool(row["veber_pass"]))
    formatted["pains_flag"] = _bool_text(bool(row["pains_flag"]))
    return formatted


def write_results(path: Path, rows: list[dict[str, object]]) -> None:
    """Write the screen table. Floats are rounded only in the file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(_csv_row(row))
