"""Figure for the drug-likeness screen.

Panel A shows each Lipinski property divided by its limit, so a bar past 1
is a violation. Panel B shows QED, coloured by the strict Lipinski call and
hatched when a PAINS alert fired. QED itself is not thresholded.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SHORT_NAMES = {
    "Epigallocatechin gallate": "EGCG",
}


def _short(name: str) -> str:
    return SHORT_NAMES.get(name, name)


def plot_screen(rows: list[dict[str, object]], path: Path) -> None:
    """Write a two-panel PNG. The file is created from ``rows`` only."""
    if not rows:
        raise ValueError("No rows to plot.")
    short = [_short(str(row["compound_name"])) for row in rows]
    series = [
        ("MW / 500", [float(row["molecular_weight"]) / 500.0 for row in rows], "#0072B2"),
        ("Crippen logP / 5", [float(row["crippen_logp"]) / 5.0 for row in rows], "#E69F00"),
        ("HBD / 5", [float(row["hbd_lipinski"]) / 5.0 for row in rows], "#009E73"),
        ("HBA / 10", [float(row["hba_lipinski"]) / 10.0 for row in rows], "#CC79A7"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.0), constrained_layout=True)
    ax = axes[0]
    indexes = list(range(len(rows)))
    width = 0.18
    for series_index, (label, values, color) in enumerate(series):
        offset = (series_index - 1.5) * width
        ax.bar(
            [index + offset for index in indexes],
            values,
            width=width,
            label=label,
            color=color,
            zorder=2,
        )
    ax.axhline(
        1.0,
        color="#222222",
        linestyle="--",
        linewidth=1.0,
        zorder=3,
        label="Lipinski limit",
    )
    ax.set_xticks(indexes)
    ax.set_xticklabels(short, rotation=30, ha="right")
    ax.set_ylabel("Property / Lipinski limit")
    ax.set_title("A   Lipinski ratios")
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, linewidth=0.4, color="#dddddd")
    ratios = [value for _, values, _ in series for value in values]
    ax.set_ylim(min(-0.55, min(ratios) - 0.15), max(2.35, max(ratios) + 0.25))
    ax.legend(frameon=False, fontsize=8, loc="upper left")

    ax2 = axes[1]
    order = sorted(range(len(rows)), key=lambda index: float(rows[index]["qed"]))
    qed_values = [float(rows[index]["qed"]) for index in order]
    colors = ["#0072B2" if rows[index]["lipinski_pass_strict"] else "#D55E00" for index in order]
    bars = ax2.barh(
        list(range(len(order))),
        qed_values,
        color=colors,
        edgecolor="#222222",
        linewidth=0.6,
        height=0.72,
        zorder=2,
    )
    for bar, index in zip(bars, order, strict=True):
        if rows[index]["pains_flag"]:
            bar.set_hatch("//")
            bar.set_hatchcolor("#ffffff")
            bar.set_hatch_linewidth(1.2)
    pains_key = Patch(facecolor="white", edgecolor="#222222", hatch="//", label="PAINS alert")
    pains_key.set_hatchcolor("#222222")
    ax2.set_yticks(list(range(len(order))))
    ax2.set_yticklabels([short[index] for index in order])
    ax2.set_xlabel("QED (RDKit, mean weights)")
    ax2.set_xlim(0, 1)
    ax2.set_title("B   QED")
    ax2.set_axisbelow(True)
    ax2.xaxis.grid(True, linewidth=0.4, color="#dddddd")
    ax2.legend(
        handles=[
            Patch(facecolor="#0072B2", edgecolor="#222222", label="No Lipinski violation"),
            Patch(facecolor="#D55E00", edgecolor="#222222", label="Lipinski violation"),
            pains_key,
        ],
        frameon=False,
        fontsize=8,
        loc="lower right",
    )

    version = rows[0]["rdkit_version"]
    fig.suptitle(
        f"Drug-likeness of {len(rows)} phytochemicals (RDKit {version})",
        fontsize=12,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)
