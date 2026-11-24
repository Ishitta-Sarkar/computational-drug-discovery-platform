"""Text summary of a drug-likeness screen."""

from dataclasses import dataclass

import rdkit

from druglikeness.filters import QED_THRESHOLD, RULE_LABELS
from druglikeness.screen import ScreenResult

_LIPINSKI_RULES = ("molecular_weight", "logp", "hbd", "hba")
_VEBER_RULES = ("tpsa", "rotatable_bonds")


@dataclass(frozen=True)
class ReferenceAgreement:
    n_ligands: int
    lipinski_flips: tuple[str, ...]
    veber_flips: tuple[str, ...]
    hbond_mismatches: tuple[str, ...]
    max_logp_gap_name: str
    max_logp_gap: float
    max_logp_rdkit: float
    max_logp_reference: float
    max_tpsa_gap_name: str
    max_tpsa_gap: float
    max_tpsa_rdkit: float
    max_tpsa_reference: float
    max_rotatable_gap: int
    max_mw_gap_name: str
    max_mw_gap: float
    max_mw_rdkit: float
    max_mw_reference: float


def reference_agreement(results: list[ScreenResult]) -> ReferenceAgreement:
    """Compare RDKit pass/fail calls with the same cutoffs on PubChem numbers.

    PubChem values never decide the screen. The comparison only shows whether
    the call for this set depends on the descriptor implementation.
    """

    lipinski_flips: list[str] = []
    veber_flips: list[str] = []
    hbond_mismatches: list[str] = []
    max_logp_name = results[0].name
    max_logp_gap = -1.0
    max_logp_rdkit = 0.0
    max_logp_reference = 0.0
    max_tpsa_name = results[0].name
    max_tpsa_gap = -1.0
    max_tpsa_rdkit = 0.0
    max_tpsa_reference = 0.0
    max_rotatable_gap = 0
    max_mw_name = results[0].name
    max_mw_gap = -1.0
    max_mw_rdkit = 0.0
    max_mw_reference = 0.0

    for result in results:
        reference = result.reference
        descriptors = result.descriptors
        ref_lipinski = _fails(
            {
                "molecular_weight": reference.molecular_weight,
                "logp": reference.xlogp,
                "hbd": float(reference.hbd),
                "hba": float(reference.hba),
            },
            {"molecular_weight": 500.0, "logp": 5.0, "hbd": 5.0, "hba": 10.0},
        )
        ref_veber = _fails(
            {
                "tpsa": reference.tpsa,
                "rotatable_bonds": float(reference.rotatable_bonds),
            },
            {"tpsa": 140.0, "rotatable_bonds": 10.0},
        )
        if set(ref_lipinski) != set(result.lipinski_violations):
            lipinski_flips.append(result.name)
        if set(ref_veber) != set(result.veber_violations):
            veber_flips.append(result.name)
        if reference.hbd != descriptors.hbd or reference.hba != descriptors.hba:
            hbond_mismatches.append(result.name)

        logp_gap = abs(descriptors.crippen_logp - reference.xlogp)
        if logp_gap > max_logp_gap:
            max_logp_gap = logp_gap
            max_logp_name = result.name
            max_logp_rdkit = descriptors.crippen_logp
            max_logp_reference = reference.xlogp
        tpsa_gap = abs(descriptors.tpsa - reference.tpsa)
        if tpsa_gap > max_tpsa_gap:
            max_tpsa_gap = tpsa_gap
            max_tpsa_name = result.name
            max_tpsa_rdkit = descriptors.tpsa
            max_tpsa_reference = reference.tpsa
        max_rotatable_gap = max(
            max_rotatable_gap,
            abs(descriptors.rotatable_bonds - reference.rotatable_bonds),
        )
        mw_gap = abs(descriptors.molecular_weight - reference.molecular_weight)
        if mw_gap > max_mw_gap:
            max_mw_gap = mw_gap
            max_mw_name = result.name
            max_mw_rdkit = descriptors.molecular_weight
            max_mw_reference = reference.molecular_weight

    return ReferenceAgreement(
        n_ligands=len(results),
        lipinski_flips=tuple(lipinski_flips),
        veber_flips=tuple(veber_flips),
        hbond_mismatches=tuple(hbond_mismatches),
        max_logp_gap_name=max_logp_name,
        max_logp_gap=max_logp_gap,
        max_logp_rdkit=max_logp_rdkit,
        max_logp_reference=max_logp_reference,
        max_tpsa_gap_name=max_tpsa_name,
        max_tpsa_gap=max_tpsa_gap,
        max_tpsa_rdkit=max_tpsa_rdkit,
        max_tpsa_reference=max_tpsa_reference,
        max_rotatable_gap=max_rotatable_gap,
        max_mw_gap_name=max_mw_name,
        max_mw_gap=max_mw_gap,
        max_mw_rdkit=max_mw_rdkit,
        max_mw_reference=max_mw_reference,
    )


def format_summary(
    results: list[ScreenResult],
    *,
    retrieved_utc: str,
    request_url: str,
) -> str:
    if not results:
        raise ValueError("Cannot summarise an empty screen.")
    agreement = reference_agreement(results)
    sections = [
        "Ligand drug-likeness screen",
        "=" * 32,
        "",
        "Question",
        "--------",
        (
            "Which of these six natural products pass Lipinski, Veber and a "
            f"QED cutoff of {QED_THRESHOLD:.1f} when every descriptor is calculated "
            "from the structure PubChem returns for its CID?"
        ),
        "",
        "Data",
        "----",
        "Source: PubChem PUG REST compound property table.",
        f"Retrieved (UTC): {retrieved_utc}",
        f"Request: {request_url}",
        f"RDKit: {rdkit.__version__}",
        (
            "The SMILES string is the only structural input. Molecular weight, "
            "Crippen logP, hydrogen-bond counts, TPSA, rotatable bonds and QED "
            "are calculated from that string. PubChem's own property numbers are "
            "stored beside them and are not used to pass or fail a ligand."
        ),
        "",
        "Rules",
        "-----",
        (
            "Lipinski, pass = zero violations: molecular weight <= 500, "
            "Crippen logP <= 5, hydrogen-bond donors <= 5, "
            "hydrogen-bond acceptors <= 10."
        ),
        "Veber, pass = zero violations: TPSA <= 140 Å², rotatable bonds <= 10.",
        (
            f"QED pass: QED >= {QED_THRESHOLD:.1f}. QED is a continuous score "
            "(Bickerton et al., 2012). The 0.5 line is a cutoff chosen for this "
            "screen, not a threshold from that paper."
        ),
        "Retained: Lipinski, Veber and the QED cutoff all pass.",
        (
            "Rotatable bonds follow Veber et al. (2002): single non-ring bonds "
            "between non-terminal heavy atoms, amide C-N bonds excluded, esters "
            "counted. RDKit's strict rotatable-bond count also drops esters and "
            "is the count inside QED. It is not the Veber column in this table."
        ),
        "",
        "Results",
        "-------",
        format_table(results),
        "",
        "Retained: " + _names(result for result in results if result.retained) + ".",
        "Not retained: " + _names(result for result in results if not result.retained) + ".",
        "",
        _lipinski_margin_note(results),
        "",
        "Rutin",
        "-----",
        _rutin_section(results),
        "",
        "Where the filters disagree",
        "--------------------------",
        _disagreement_section(results),
        "",
        "PubChem property check",
        "----------------------",
        _agreement_section(agreement),
        "",
        "Scope",
        "-----",
        _scope_section(),
    ]
    return "\n".join(sections).rstrip()


def format_table(results: list[ScreenResult]) -> str:
    headers = [
        "Ligand",
        "MW",
        "logP",
        "HBD",
        "HBA",
        "TPSA",
        "RotB",
        "QED",
        "Lipinski",
        "Veber",
        "QED>=0.5",
        "Retained",
    ]
    rows = [headers]
    for result in results:
        descriptors = result.descriptors
        rows.append(
            [
                result.name,
                f"{descriptors.molecular_weight:.2f}",
                f"{descriptors.crippen_logp:.2f}",
                str(descriptors.hbd),
                str(descriptors.hba),
                f"{descriptors.tpsa:.2f}",
                str(descriptors.rotatable_bonds),
                f"{descriptors.qed:.3f}",
                "pass" if result.lipinski_pass else "fail",
                "pass" if result.veber_pass else "fail",
                "pass" if result.qed_pass else "fail",
                "yes" if result.retained else "no",
            ]
        )
    widths = [max(len(row[index]) for row in rows) for index in range(len(headers))]
    lines = []
    for row_index, row in enumerate(rows):
        line = "  ".join(value.ljust(widths[index]) for index, value in enumerate(row))
        lines.append(line.rstrip())
        if row_index == 0:
            lines.append("  ".join("-" * width for width in widths))
    return "\n".join(lines)


def _rutin_section(results: list[ScreenResult]) -> str:
    rutin = _by_name(results, "Rutin")
    quercetin = _by_name(results, "Quercetin")
    if rutin is None or quercetin is None:
        return "Rutin or quercetin is absent from this screen, so the comparison was skipped."
    rutin_d = rutin.descriptors
    quercetin_d = quercetin.descriptors
    rutin_rules = _join_labels(rutin.lipinski_violations)
    return "\n".join(
        [
            (
                "Rutin is quercetin-3-O-rutinoside: the flavonol quercetin with a "
                "rutinose disaccharide on the 3-hydroxy group. The two rows are "
                "the same core with and without that sugar."
            ),
            (
                f"Quercetin passes Lipinski (MW {quercetin_d.molecular_weight:.2f}, "
                f"Crippen logP {quercetin_d.crippen_logp:.2f}, "
                f"HBD {quercetin_d.hbd}, HBA {quercetin_d.hba}) and Veber "
                f"(TPSA {quercetin_d.tpsa:.2f}, rotatable bonds "
                f"{quercetin_d.rotatable_bonds})."
            ),
            (
                f"Rutin fails Lipinski on {rutin_rules}. "
                f"MW is {rutin_d.molecular_weight:.2f} (cutoff 500), "
                f"hydrogen-bond donors are {rutin_d.hbd} (cutoff 5) and "
                f"acceptors are {rutin_d.hba} (cutoff 10). "
                f"Crippen logP is {rutin_d.crippen_logp:.2f}, which is inside the "
                "Lipinski logP limit. The failures are size and hydrogen bonding, "
                "not lipophilicity."
            ),
            (
                f"Veber fails on TPSA ({rutin_d.tpsa:.2f} against a cutoff of 140). "
                f"Rotatable bonds are {rutin_d.rotatable_bonds}, within the limit of 10. "
                f"QED falls from {quercetin_d.qed:.3f} for quercetin to "
                f"{rutin_d.qed:.3f} for rutin."
            ),
            (
                "The disaccharide is what crosses the cutoffs. These rules compare "
                "the intact structure with property ranges associated with oral drugs. "
                "They do not measure whether rutin binds a protein, and they do not "
                "say whether the sugar is cleaved before absorption. A glycoside can "
                "be a reasonable natural-product sample and a poor match to Lipinski's "
                "oral ranges at the same time."
            ),
            (
                "An earlier version of this repository sorted the same six names by "
                "one docking_score column and printed rutin first. The column had no "
                "target, and the summary never said that rutin fails Lipinski. That "
                "column has been removed. A docking score answers a different question "
                "from these filters, and only after a documented ligand-target calculation."
            ),
        ]
    )


def _disagreement_section(results: list[ScreenResult]) -> str:
    paragraphs = []
    egcg = _by_name(results, "Epigallocatechin gallate")
    if egcg is not None and not egcg.lipinski_pass:
        descriptors = egcg.descriptors
        paragraphs.append(
            "Epigallocatechin gallate is the other Lipinski failure, and it is "
            f"not a glycoside. RDKit counts {descriptors.hbd} donors and "
            f"{descriptors.hba} acceptors, TPSA is {descriptors.tpsa:.2f}, and "
            f"QED is {descriptors.qed:.3f}. The gallate and the polyphenol pattern "
            "are enough to leave the oral ranges. Passing is not the default for "
            "a well-known natural product in this set."
        )
    quercetin = _by_name(results, "Quercetin")
    luteolin = _by_name(results, "Luteolin")
    curcumin = _by_name(results, "Curcumin")
    if (
        quercetin is not None
        and luteolin is not None
        and curcumin is not None
        and quercetin.lipinski_pass
        and quercetin.veber_pass
        and not quercetin.qed_pass
    ):
        paragraphs.append(
            "Quercetin is the awkward middle case. It passes Lipinski and Veber, "
            f"but QED is {quercetin.descriptors.qed:.3f}, below the "
            f"{QED_THRESHOLD:.1f} cutoff used here, so it is not retained. "
            f"Luteolin, which lacks the 3-hydroxy group, has QED "
            f"{luteolin.descriptors.qed:.3f} and is retained. Curcumin has QED "
            f"{curcumin.descriptors.qed:.3f}. Those three scores sit close to a "
            "line this workflow chose. They should not be read as a sharp "
            "separation into drugs and non-drugs."
        )
    if not paragraphs:
        return "No extra disagreement to report for this set."
    return "\n".join(paragraphs)


def _lipinski_margin_note(results: list[ScreenResult]) -> str:
    single = [result.name for result in results if len(result.lipinski_violations) == 1]
    if single:
        return (
            "Ligands with exactly one Lipinski violation: "
            + ", ".join(single)
            + ". The rule of five is often applied with one violation allowed. "
            "This workflow still requires zero violations."
        )
    return (
        "No ligand in this set has exactly one Lipinski violation. "
        "Allowing one violation, which is how the rule of five is often applied, "
        "would not change who passes Lipinski here. The workflow still requires "
        "zero violations, and the violation list shows how large each failure is."
    )


def _agreement_section(agreement: ReferenceAgreement) -> str:
    if agreement.lipinski_flips or agreement.veber_flips:
        lipinski = ", ".join(agreement.lipinski_flips) or "none"
        veber = ", ".join(agreement.veber_flips) or "none"
        decision = (
            "Applying the same numeric cutoffs to PubChem's values changes the "
            f"Lipinski call for: {lipinski}. It changes the Veber call for: {veber}."
        )
    else:
        decision = (
            f"Applying the same numeric cutoffs to PubChem's XLogP, TPSA, "
            f"hydrogen-bond counts and rotatable-bond counts gives the same "
            f"Lipinski call and the same Veber call for all {agreement.n_ligands} ligands."
        )
    if agreement.hbond_mismatches:
        hbond = "Hydrogen-bond counts differ for: " + ", ".join(agreement.hbond_mismatches) + "."
    else:
        hbond = "Hydrogen-bond donor and acceptor counts match PubChem for every ligand."
    if agreement.max_rotatable_gap == 0:
        rotatable = "Rotatable-bond counts match PubChem for every ligand."
    else:
        rotatable = (
            f"Rotatable-bond counts differ from PubChem by at most {agreement.max_rotatable_gap}."
        )
    return " ".join(
        [
            decision,
            hbond,
            (
                f"The largest Crippen-versus-XLogP gap is {agreement.max_logp_gap:.2f} "
                f"for {agreement.max_logp_gap_name} "
                f"({agreement.max_logp_rdkit:.2f} versus {agreement.max_logp_reference:.2f})."
            ),
            (
                f"The largest TPSA gap is {agreement.max_tpsa_gap:.2f} Å² for "
                f"{agreement.max_tpsa_gap_name} "
                f"({agreement.max_tpsa_rdkit:.2f} versus {agreement.max_tpsa_reference:.2f})."
            ),
            rotatable,
            (
                "The largest average-mass gap is "
                f"{agreement.max_mw_gap:.3f} Da for {agreement.max_mw_gap_name} "
                f"({agreement.max_mw_rdkit:.3f} versus {agreement.max_mw_reference:.3f})."
            ),
            "None of those reference numbers is an input to the decision.",
        ]
    )


def _scope_section() -> str:
    return "\n".join(
        [
            (
                "This is a ligand drug-likeness screen of a fixed six-compound "
                "learning set. It does not dock, it does not rank binding, and it "
                "does not claim activity against ACE2, PPARG, APP, SNCA or any "
                "other protein. A docking score is defined for one ligand-target "
                "pair from a documented calculation. No such calculation is in "
                "this repository."
            ),
            (
                "Lipinski's 1997 analysis set aside antibiotics, antifungals, "
                "vitamins and cardiac glycosides, many of which are taken up by "
                "transporters. A glycosylated flavonoid falling outside the cutoffs "
                "is consistent with that scope. The filters estimate resemblance to "
                "oral drugs. They are not a reason, by themselves, to discard a "
                "natural product from biological study."
            ),
        ]
    )


def _by_name(results: list[ScreenResult], name: str) -> ScreenResult | None:
    folded = name.casefold()
    for result in results:
        if result.name.casefold() == folded:
            return result
    return None


def _names(results) -> str:
    names = [result.name for result in results]
    if not names:
        return "none"
    return ", ".join(names)


def _fails(values: dict[str, float], limits: dict[str, float]) -> tuple[str, ...]:
    return tuple(rule for rule, limit in limits.items() if values[rule] > limit)


def _join_labels(rules: tuple[str, ...]) -> str:
    labels = [
        RULE_LABELS[rule] for rule in rules if rule in _LIPINSKI_RULES or rule in _VEBER_RULES
    ]
    if not labels:
        return "no listed rule"
    if len(labels) == 1:
        return labels[0]
    if len(labels) == 2:
        return f"{labels[0]} and {labels[1]}"
    return ", ".join(labels[:-1]) + " and " + labels[-1]
