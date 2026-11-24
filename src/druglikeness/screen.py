"""Apply the drug-likeness filters to retrieved structures."""

import csv
from dataclasses import dataclass
from pathlib import Path

from druglikeness.descriptors import DescriptorSet, compute_descriptors
from druglikeness.filters import lipinski_violations, passes_qed, veber_violations
from druglikeness.pubchem import CompoundRecord, ReferenceProperties


@dataclass(frozen=True)
class ScreenResult:
    name: str
    pubchem_cid: int
    smiles: str
    molecular_formula: str
    iupac_name: str
    descriptors: DescriptorSet
    reference: ReferenceProperties
    lipinski_violations: tuple[str, ...]
    veber_violations: tuple[str, ...]

    @property
    def lipinski_pass(self) -> bool:
        return not self.lipinski_violations

    @property
    def veber_pass(self) -> bool:
        return not self.veber_violations

    @property
    def qed_pass(self) -> bool:
        return passes_qed(self.descriptors)

    @property
    def retained(self) -> bool:
        return self.lipinski_pass and self.veber_pass and self.qed_pass


def screen_records(
    records: list[CompoundRecord] | tuple[CompoundRecord, ...],
) -> list[ScreenResult]:
    results: list[ScreenResult] = []
    for record in records:
        descriptors = compute_descriptors(record.smiles, record.molecular_formula)
        results.append(
            ScreenResult(
                name=record.name,
                pubchem_cid=record.pubchem_cid,
                smiles=record.smiles,
                molecular_formula=record.molecular_formula,
                iupac_name=record.iupac_name,
                descriptors=descriptors,
                reference=record.reference,
                lipinski_violations=lipinski_violations(descriptors),
                veber_violations=veber_violations(descriptors),
            )
        )
    return results


TABLE_FIELDS = (
    "compound_name",
    "pubchem_cid",
    "molecular_formula",
    "smiles",
    "molecular_weight",
    "crippen_logp",
    "hbd",
    "hba",
    "tpsa",
    "rotatable_bonds",
    "qed",
    "lipinski_violations",
    "veber_violations",
    "lipinski_pass",
    "veber_pass",
    "qed_pass",
    "retained",
    "pubchem_molecular_weight",
    "pubchem_xlogp",
    "pubchem_tpsa",
    "pubchem_hbd",
    "pubchem_hba",
    "pubchem_rotatable_bonds",
)


def write_table(path: Path, results: list[ScreenResult]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TABLE_FIELDS)
        writer.writeheader()
        for result in results:
            descriptors = result.descriptors
            reference = result.reference
            writer.writerow(
                {
                    "compound_name": result.name,
                    "pubchem_cid": result.pubchem_cid,
                    "molecular_formula": result.molecular_formula,
                    "smiles": result.smiles,
                    "molecular_weight": f"{descriptors.molecular_weight:.10f}",
                    "crippen_logp": f"{descriptors.crippen_logp:.10f}",
                    "hbd": descriptors.hbd,
                    "hba": descriptors.hba,
                    "tpsa": f"{descriptors.tpsa:.10f}",
                    "rotatable_bonds": descriptors.rotatable_bonds,
                    "qed": f"{descriptors.qed:.10f}",
                    "lipinski_violations": ";".join(result.lipinski_violations),
                    "veber_violations": ";".join(result.veber_violations),
                    "lipinski_pass": _bool(result.lipinski_pass),
                    "veber_pass": _bool(result.veber_pass),
                    "qed_pass": _bool(result.qed_pass),
                    "retained": _bool(result.retained),
                    "pubchem_molecular_weight": reference.molecular_weight,
                    "pubchem_xlogp": reference.xlogp,
                    "pubchem_tpsa": reference.tpsa,
                    "pubchem_hbd": reference.hbd,
                    "pubchem_hba": reference.hba,
                    "pubchem_rotatable_bonds": reference.rotatable_bonds,
                }
            )


def _bool(value: bool) -> str:
    return "true" if value else "false"
