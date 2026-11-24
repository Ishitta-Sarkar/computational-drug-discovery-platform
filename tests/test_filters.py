from druglikeness.descriptors import DescriptorSet
from druglikeness.filters import (
    QED_THRESHOLD,
    lipinski_violations,
    passes_qed,
    veber_violations,
)


def _descriptors(**overrides) -> DescriptorSet:
    values = {
        "molecular_weight": 300.0,
        "crippen_logp": 2.0,
        "hbd": 2,
        "hba": 4,
        "tpsa": 80.0,
        "rotatable_bonds": 3,
        "qed": 0.6,
    }
    values.update(overrides)
    return DescriptorSet(**values)


def test_each_lipinski_rule_can_fail_on_its_own():
    assert lipinski_violations(_descriptors(molecular_weight=500.0)) == ()
    assert lipinski_violations(_descriptors(molecular_weight=500.1)) == ("molecular_weight",)
    assert lipinski_violations(_descriptors(crippen_logp=5.0)) == ()
    assert lipinski_violations(_descriptors(crippen_logp=5.1)) == ("logp",)
    assert lipinski_violations(_descriptors(hbd=5)) == ()
    assert lipinski_violations(_descriptors(hbd=6)) == ("hbd",)
    assert lipinski_violations(_descriptors(hba=10)) == ()
    assert lipinski_violations(_descriptors(hba=11)) == ("hba",)


def test_negative_logp_is_inside_the_lipinski_limit():
    assert "logp" not in lipinski_violations(_descriptors(crippen_logp=-1.69))


def test_veber_boundaries():
    assert veber_violations(_descriptors(tpsa=140.0, rotatable_bonds=10)) == ()
    assert veber_violations(_descriptors(tpsa=140.1)) == ("tpsa",)
    assert veber_violations(_descriptors(rotatable_bonds=11)) == ("rotatable_bonds",)


def test_qed_cutoff_is_inclusive():
    assert passes_qed(_descriptors(qed=QED_THRESHOLD))
    assert not passes_qed(_descriptors(qed=QED_THRESHOLD - 0.001))
