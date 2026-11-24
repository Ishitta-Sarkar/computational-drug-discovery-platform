"""Results figure for the drug-likeness screen."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from druglikeness.filters import LIPINSKI_LIMITS, QED_THRESHOLD, VEBER_LIMITS
from druglikeness.screen import ScreenResult

WITHIN = "#D6E6F5"
ABOVE = "#F6D6C9"
RETAINED = "#2C6EAA"
EXCLUDED = "#C65A2E"
TEXT = "#1A1A1A"

_COLUMNS = (
    ("molecular_weight", "MW\n≤ 500", "{:.1f}"),
    ("logp", "logP\n≤ 5", "{:.2f}"),
    ("hbd", "HBD\n≤ 5", "{:.0f}"),
    ("hba", "HBA\n≤ 10", "{:.0f}"),
    ("tpsa", "TPSA\n≤ 140", "{:.1f}"),
    ("rotatable_bonds", "RotB\n≤ 10", "{:.0f}"),
)


def plot_screen(results: list[ScreenResult], path: Path) -> None:
    """Write a two-panel figure: rule checks, then QED.

    Ligands are ordered by QED, highest at the top. Cell color marks whether
    a Lipinski or Veber value is within its cutoff. Bar color marks the
    combined decision, so a ligand can be entirely within the structural
    cutoffs and still be drawn as not retained when QED is below the cutoff.
    """

    if not results:
        raise ValueError("Cannot plot an empty screen.")
    ordered = sorted(results, key=lambda result: (-result.descriptors.qed, result.name))
    names = [result.name for result in ordered]
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(11.4, 5.8),
        gridspec_kw={"width_ratios": [1.45, 1.0]},
        constrained_layout=True,
    )
    _draw_rules(axes[0], ordered, names)
    _draw_qed(axes[1], ordered, names)
    figure.suptitle("Ligand drug-likeness screen", fontsize=13, color=TEXT)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180, facecolor="white")
    plt.close(figure)


def _draw_rules(axis, results: list[ScreenResult], names: list[str]) -> None:
    status = []
    annotations = []
    for result in results:
        values = _rule_values(result)
        row_status = []
        row_text = []
        for key, _label, template in _COLUMNS:
            value = values[key]
            limit = _limit_for(key)
            row_status.append(0 if value <= limit else 1)
            row_text.append(template.format(value))
        status.append(row_status)
        annotations.append(row_text)

    axis.imshow(
        np.array(status),
        cmap=ListedColormap([WITHIN, ABOVE]),
        vmin=0,
        vmax=1,
        aspect="auto",
    )
    axis.set_xticks(range(len(_COLUMNS)))
    axis.set_xticklabels([column[1] for column in _COLUMNS], fontsize=8)
    axis.set_yticks(range(len(names)))
    axis.set_yticklabels(names, fontsize=9)
    for row_index, row_text in enumerate(annotations):
        for column_index, text in enumerate(row_text):
            axis.text(
                column_index,
                row_index,
                text,
                ha="center",
                va="center",
                fontsize=8,
                color=TEXT,
            )
    axis.set_xticks(np.arange(len(_COLUMNS) + 1) - 0.5, minor=True)
    axis.set_yticks(np.arange(len(names) + 1) - 0.5, minor=True)
    axis.grid(which="minor", color="white", linestyle="-", linewidth=2)
    axis.tick_params(which="minor", bottom=False, left=False)
    axis.tick_params(which="major", length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.set_title("A   Lipinski and Veber", loc="left", fontsize=11, color=TEXT, pad=8)
    axis.legend(
        handles=[
            Patch(facecolor=WITHIN, edgecolor="#7A93A8", label="Within cutoff"),
            Patch(facecolor=ABOVE, edgecolor="#C4A090", label="Above cutoff"),
        ],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=2,
        fontsize=8,
    )


def _draw_qed(axis, results: list[ScreenResult], names: list[str]) -> None:
    positions = np.arange(len(results))
    values = [result.descriptors.qed for result in results]
    colors = [RETAINED if result.retained else EXCLUDED for result in results]
    axis.barh(positions, values, color=colors, height=0.62, zorder=2)
    axis.axvline(QED_THRESHOLD, color="#222222", linestyle="--", linewidth=0.9, zorder=1)
    for position, value in zip(positions, values, strict=True):
        axis.text(
            1.02,
            position,
            f"{value:.3f}",
            va="center",
            ha="left",
            fontsize=8,
            color=TEXT,
        )
    axis.set_yticks(positions)
    axis.set_yticklabels(names, fontsize=9)
    axis.invert_yaxis()
    axis.set_xlim(0, 1.28)
    axis.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    axis.set_xlabel("QED")
    axis.set_title("B   QED", loc="left", fontsize=11, color=TEXT, pad=8)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.legend(
        handles=[
            Patch(facecolor=RETAINED, label="Retained"),
            Patch(facecolor=EXCLUDED, label="Not retained"),
            plt.Line2D([0], [0], color="#222222", linestyle="--", label="Cutoff 0.5"),
        ],
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=3,
        fontsize=8,
    )


def _rule_values(result: ScreenResult) -> dict[str, float]:
    descriptors = result.descriptors
    return {
        "molecular_weight": descriptors.molecular_weight,
        "logp": descriptors.crippen_logp,
        "hbd": float(descriptors.hbd),
        "hba": float(descriptors.hba),
        "tpsa": descriptors.tpsa,
        "rotatable_bonds": float(descriptors.rotatable_bonds),
    }


def _limit_for(rule: str) -> float:
    if rule in LIPINSKI_LIMITS:
        return LIPINSKI_LIMITS[rule]
    return VEBER_LIMITS[rule]
