"""Lipinski, Veber and QED decision rules.

Cutoffs
-------
Lipinski (zero violations required by this workflow):
    molecular weight <= 500
    Crippen logP <= 5
    hydrogen-bond donors <= 5
    hydrogen-bond acceptors <= 10

Veber:
    TPSA <= 140
    rotatable bonds <= 10

QED:
    >= 0.5

The QED cutoff is a choice of this workflow. Bickerton et al. published QED
as a continuous desirability score and did not define a universal pass line.
Lipinski's rule is often applied with one violation allowed. This workflow
requires zero Lipinski violations and reports the violation list so a miss
by one rule is visible separately from a miss by three.
"""

from druglikeness.descriptors import DescriptorSet

LIPINSKI_LIMITS: dict[str, float] = {
    "molecular_weight": 500.0,
    "logp": 5.0,
    "hbd": 5.0,
    "hba": 10.0,
}

VEBER_LIMITS: dict[str, float] = {
    "tpsa": 140.0,
    "rotatable_bonds": 10.0,
}

QED_THRESHOLD = 0.5

RULE_LABELS: dict[str, str] = {
    "molecular_weight": "molecular weight",
    "logp": "Crippen logP",
    "hbd": "hydrogen-bond donors",
    "hba": "hydrogen-bond acceptors",
    "tpsa": "TPSA",
    "rotatable_bonds": "rotatable bonds",
}


def lipinski_values(descriptors: DescriptorSet) -> dict[str, float]:
    return {
        "molecular_weight": descriptors.molecular_weight,
        "logp": descriptors.crippen_logp,
        "hbd": float(descriptors.hbd),
        "hba": float(descriptors.hba),
    }


def veber_values(descriptors: DescriptorSet) -> dict[str, float]:
    return {
        "tpsa": descriptors.tpsa,
        "rotatable_bonds": float(descriptors.rotatable_bonds),
    }


def upper_bound_violations(values: dict[str, float], limits: dict[str, float]) -> tuple[str, ...]:
    violations = []
    for rule, limit in limits.items():
        if values[rule] > limit:
            violations.append(rule)
    return tuple(violations)


def lipinski_violations(descriptors: DescriptorSet) -> tuple[str, ...]:
    return upper_bound_violations(lipinski_values(descriptors), LIPINSKI_LIMITS)


def veber_violations(descriptors: DescriptorSet) -> tuple[str, ...]:
    return upper_bound_violations(veber_values(descriptors), VEBER_LIMITS)


def passes_qed(descriptors: DescriptorSet, threshold: float = QED_THRESHOLD) -> bool:
    return descriptors.qed >= threshold
