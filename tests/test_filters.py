"""Rule thresholds, independent of the ligand cache."""

from druglikeness.descriptors import MolecularDescriptors
from druglikeness.filters import (
    LIPINSKI_MAX_HBA,
    LIPINSKI_MAX_HBD,
    LIPINSKI_MAX_LOGP,
    LIPINSKI_MAX_MOLECULAR_WEIGHT,
    VEBER_MAX_ROTATABLE_BONDS,
    VEBER_MAX_TPSA,
    lipinski_alert,
    lipinski_assessment,
    veber_assessment,
)


def _desc(**overrides: float) -> MolecularDescriptors:
    values = {
        "molecular_weight": 300.0,
        "crippen_logp": 2.0,
        "hbd_lipinski": 2,
        "hba_lipinski": 4,
        "tpsa": 80.0,
        "rotatable_bonds": 3,
        "qed": 0.5,
        "inchikey": "TEST",
    }
    values.update(overrides)
    return MolecularDescriptors(**values)


def test_values_on_the_lipinski_limits_pass() -> None:
    assessment = lipinski_assessment(
        _desc(
            molecular_weight=LIPINSKI_MAX_MOLECULAR_WEIGHT,
            crippen_logp=LIPINSKI_MAX_LOGP,
            hbd_lipinski=LIPINSKI_MAX_HBD,
            hba_lipinski=LIPINSKI_MAX_HBA,
        )
    )
    assert assessment.passed
    assert assessment.count == 0
    assert lipinski_alert(assessment) is False


def test_one_lipinski_violation_fails_strict_screen_and_does_not_alert() -> None:
    assessment = lipinski_assessment(_desc(molecular_weight=500.01))
    assert assessment.violations == ("molecular_weight",)
    assert not assessment.passed
    assert lipinski_alert(assessment) is False


def test_two_lipinski_violations_raise_the_original_alert() -> None:
    assessment = lipinski_assessment(_desc(hbd_lipinski=6, hba_lipinski=11))
    assert assessment.violations == ("hbd_lipinski", "hba_lipinski")
    assert lipinski_alert(assessment) is True
    assert assessment.labels() == "H-bond donors; H-bond acceptors"


def test_veber_limits_are_inclusive() -> None:
    assessment = veber_assessment(
        _desc(rotatable_bonds=VEBER_MAX_ROTATABLE_BONDS, tpsa=VEBER_MAX_TPSA)
    )
    assert assessment.passed


def test_veber_fails_on_tpsa_alone() -> None:
    assessment = veber_assessment(_desc(tpsa=VEBER_MAX_TPSA + 0.01, rotatable_bonds=1))
    assert assessment.violations == ("tpsa",)
    assert not assessment.passed
