"""Plain-text summary of a screen.

The rutin paragraph is built from the computed violation lists. It is
included only when both rutin and quercetin are in the table, because the
comparison is between that glycoside and its aglycone.
"""

from __future__ import annotations

from druglikeness.filters import VIOLATION_LABELS

_SHORT_LABELS = {
    "molecular_weight": "MW",
    "crippen_logp": "logP",
    "hbd_lipinski": "HBD",
    "hba_lipinski": "HBA",
    "rotatable_bonds": "RB",
    "tpsa": "TPSA",
}


def _labels(violation_field: str, short: bool = False) -> str:
    if not violation_field:
        return "pass"
    keys = [key for key in violation_field.split(";") if key]
    names = _SHORT_LABELS if short else VIOLATION_LABELS
    return "; ".join(names[key] for key in keys)


def _row_by_name(rows: list[dict[str, object]], name: str) -> dict[str, object] | None:
    for row in rows:
        if str(row["compound_name"]).casefold() == name.casefold():
            return row
    return None


def format_report(rows: list[dict[str, object]], retrieved_on: str) -> str:
    """Return the human-readable report printed by ``main.py``."""
    if not rows:
        raise ValueError("No rows to report.")
    rdkit_version = rows[0]["rdkit_version"]
    headers = (
        "Compound",
        "MW",
        "logP",
        "HBD",
        "HBA",
        "TPSA",
        "RB",
        "QED",
        "Lipinski",
        "Veber",
        "PAINS",
    )
    table: list[tuple[str, ...]] = [headers]
    for row in rows:
        pains = str(row["pains_alerts"]) if row["pains_alerts"] else "none"
        table.append(
            (
                str(row["compound_name"]),
                f"{float(row['molecular_weight']):.2f}",
                f"{float(row['crippen_logp']):.2f}",
                str(row["hbd_lipinski"]),
                str(row["hba_lipinski"]),
                f"{float(row['tpsa']):.2f}",
                str(row["rotatable_bonds"]),
                f"{float(row['qed']):.3f}",
                _labels(str(row["lipinski_violations"]), short=True),
                _labels(str(row["veber_violations"]), short=True),
                pains,
            )
        )
    widths = [max(len(record[index]) for record in table) for index in range(len(headers))]
    lines = [
        "Ligand drug-likeness screen",
        f"RDKit {rdkit_version}",
        f"Compounds: {len(rows)}",
        f"PubChem cache retrieved: {retrieved_on}",
        "",
    ]
    for record in table:
        lines.append("  ".join(value.ljust(widths[i]) for i, value in enumerate(record)).rstrip())
    lines.extend(
        [
            "",
            "Lipinski pass means zero violations "
            "(MW <= 500, Crippen logP <= 5, HBD <= 5, HBA <= 10).",
            "A Lipinski alert, following the 1997 wording, is two or more violations.",
            "Veber pass means rotatable bonds <= 10 and TPSA <= 140 Å².",
            "QED is a score from 0 to 1. It is not turned into a pass or fail.",
            "PAINS alerts are substructure flags, not activity calls.",
            "",
        ]
    )
    lines.extend(_rutin_lines(rows))
    lines.append(
        "No docking score is reported. A docking score needs a ligand, "
        "a target structure and a docking run, and this repository does not run docking."
    )
    return "\n".join(lines)


def _rutin_lines(rows: list[dict[str, object]]) -> list[str]:
    rutin = _row_by_name(rows, "Rutin")
    quercetin = _row_by_name(rows, "Quercetin")
    if rutin is None or quercetin is None:
        return []
    rutin_lip = _labels(str(rutin["lipinski_violations"]))
    rutin_veber = _labels(str(rutin["veber_violations"]))
    quercetin_lip = _labels(str(quercetin["lipinski_violations"]))
    quercetin_veber = _labels(str(quercetin["veber_violations"]))
    return [
        "Rutin compared with quercetin, its aglycone:",
        f"  Rutin Lipinski: {rutin_lip}. Rutin Veber: {rutin_veber}.",
        f"  Quercetin Lipinski: {quercetin_lip}. Quercetin Veber: {quercetin_veber}.",
        "  The earlier table in this repository ranked rutin first with a docking "
        "score of -9.0 kcal/mol and no target. That score was not produced by a "
        "run stored here, and it is not used.",
        "",
    ]
