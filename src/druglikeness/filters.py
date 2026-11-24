"""Lipinski, Veber and PAINS calls on RDKit descriptors.

Lipinski et al. (1997) described poor absorption or permeation as more likely
when two or more of these hold: molecular weight over 500, logP over 5, more
than 5 hydrogen-bond donors, more than 10 hydrogen-bond acceptors. A strict
screen, the one used by the first version of this repository, keeps only
molecules with zero violations. Both calls are reported.

Veber et al. (2002) associated oral bioavailability in their rat data with
10 or fewer rotatable bonds and a polar surface area of 140 Å² or less.
Both conditions are required for a Veber pass.

Cutoffs use '>' so a value sitting on the published limit still passes.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

from druglikeness.descriptors import MolecularDescriptors

LIPINSKI_MAX_MOLECULAR_WEIGHT = 500.0
LIPINSKI_MAX_LOGP = 5.0
LIPINSKI_MAX_HBD = 5
LIPINSKI_MAX_HBA = 10

VEBER_MAX_ROTATABLE_BONDS = 10
VEBER_MAX_TPSA = 140.0

VIOLATION_LABELS = {
    "molecular_weight": "molecular weight",
    "crippen_logp": "Crippen logP",
    "hbd_lipinski": "H-bond donors",
    "hba_lipinski": "H-bond acceptors",
    "rotatable_bonds": "rotatable bonds",
    "tpsa": "TPSA",
}


@dataclass(frozen=True)
class RuleAssessment:
    """Which named criteria a molecule missed."""

    violations: tuple[str, ...]

    @property
    def count(self) -> int:
        return len(self.violations)

    @property
    def passed(self) -> bool:
        return self.count == 0

    def labels(self) -> str:
        if not self.violations:
            return "pass"
        return "; ".join(VIOLATION_LABELS[key] for key in self.violations)


def lipinski_assessment(desc: MolecularDescriptors) -> RuleAssessment:
    """Return Lipinski violations. An empty tuple is a strict pass."""
    violations: list[str] = []
    if desc.molecular_weight > LIPINSKI_MAX_MOLECULAR_WEIGHT:
        violations.append("molecular_weight")
    if desc.crippen_logp > LIPINSKI_MAX_LOGP:
        violations.append("crippen_logp")
    if desc.hbd_lipinski > LIPINSKI_MAX_HBD:
        violations.append("hbd_lipinski")
    if desc.hba_lipinski > LIPINSKI_MAX_HBA:
        violations.append("hba_lipinski")
    return RuleAssessment(tuple(violations))


def lipinski_alert(assessment: RuleAssessment) -> bool:
    """True when two or more Lipinski parameters are out of range.

    This matches the original wording more closely than a zero-violation
    screen. One violation does not raise the alert.
    """
    return assessment.count >= 2


def veber_assessment(desc: MolecularDescriptors) -> RuleAssessment:
    """Return Veber violations. An empty tuple is a pass."""
    violations: list[str] = []
    if desc.rotatable_bonds > VEBER_MAX_ROTATABLE_BONDS:
        violations.append("rotatable_bonds")
    if desc.tpsa > VEBER_MAX_TPSA:
        violations.append("tpsa")
    return RuleAssessment(tuple(violations))


def _pains_catalog() -> FilterCatalog.FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog.FilterCatalog(params)


def pains_descriptions(mol: Chem.Mol) -> list[str]:
    """RDKit PAINS filter descriptions matched by this molecule.

    The strings are the catalog names, for example ``catechol_A(92)``.
    The letter inside the name comes from the Baell and Holloway filter
    label. It is not a separate assay result. An empty list means the
    combined PAINS catalog did not match.
    """
    found: list[str] = []
    for entry in _pains_catalog().GetMatches(mol):
        description = entry.GetDescription()
        if description not in found:
            found.append(description)
    return found
