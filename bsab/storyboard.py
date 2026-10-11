"""Figures built to *explain* the project to a non-specialist audience.

`report.py` makes figures that evidence the results. These make the idea
legible to someone who has never heard of a bispecific antibody: the pipeline
as a funnel, why same-cell co-expression matters, and one worked target
(BCMA in myeloma) that shows the whole argument end to end.

Chart rules are the project's usual ones: one axis per chart, emphasis colour
rather than categorical where a single series is the point, direct labels
instead of a value on every mark, and a single-hue sequential ramp for
magnitude.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import config as C
from . import data_hpa as D
from .report import (AQUA, BLUE, CRITICAL, DEEMPH, GOOD, GRID, INK, INK_2,
                     MUTED, ORANGE, SURFACE, _style)


def fig_pipeline(universe_n: int, pairs_n: int, malignant_cells: int,
                 path=None) -> None:
    """The five stages as labelled boxes.

    Deliberately not a bar chart: gene counts and pair counts are different
    units, so putting them on one axis would imply a comparison that is not
    meaningful.
    """
    path = path or (C.FIGURES / "story_pipeline.png")
    stages = [
        ("20,151", "human genes", "everything measured by\nthe Human Protein Atlas", DEEMPH),
        (f"{universe_n:,}", "an antibody can reach", "on the cell surface, so a drug\ncan physically bind them", BLUE),
        (f"{pairs_n:,}", "candidate pairs", "every sensible two-target\ncombination", BLUE),
        (f"{malignant_cells:,}", "tumour cells checked", "are both targets on the\nSAME cancer cell?", ORANGE),
        ("64", "known drugs as a test", "does the method rediscover\nbispecifics that already work?", GOOD),
    ]
    fig, ax = plt.subplots(figsize=(13.5, 3.4))
    ax.set_xlim(0, len(stages)); ax.set_ylim(0, 1); ax.axis("off")

    for i, (big, label, sub, colour) in enumerate(stages):
        x = i + 0.5
        box = mpatches.FancyBboxPatch(
            (i + 0.06, 0.18), 0.88, 0.66, boxstyle="round,pad=0.015,rounding_size=0.03",
            linewidth=1.6, edgecolor=colour, facecolor=SURFACE, zorder=2)
        ax.add_patch(box)
        ax.text(x, 0.70, big, ha="center", va="center", fontsize=19,
                fontweight="bold", color=colour, zorder=3)
        ax.text(x, 0.555, label, ha="center", va="center", fontsize=10.5,
                color=INK, zorder=3)
        ax.text(x, 0.36, sub, ha="center", va="center", fontsize=8.4,
                color=MUTED, linespacing=1.45, zorder=3)
        if i < len(stages) - 1:
            ax.annotate("", xy=(i + 1.05, 0.51), xytext=(i + 0.96, 0.51),
                        arrowprops=dict(arrowstyle="-|>", color=INK_2, lw=1.4))
    ax.set_title("How the pipeline narrows the problem", loc="left",
                 fontsize=13, fontweight="bold", color=INK, pad=12)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def fig_same_cell_concept(path=None) -> None:
    """Why averaged data is not enough - the core idea, drawn.

    A bispecific that needs both antigens has to find them on one cell. Two
    tumours can look identical on averaged expression and behave completely
    differently, and only single-cell data separates them.
    """
    path = path or (C.FIGURES / "story_same_cell.png")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))

    panels = [
        ("WORKS: both targets on the same cell", GOOD, True),
        ("FAILS: targets on different cells", CRITICAL, False),
    ]
    centres = [(2.6, 7.4), (7.2, 7.8), (4.9, 4.9), (2.4, 2.4), (7.4, 2.6)]
    radius = 1.5
    for ax, (title, colour, same) in zip(axes, panels):
        ax.set_xlim(0, 10); ax.set_ylim(-0.6, 10)
        ax.set_aspect("equal")  # otherwise the cells render as ellipses
        ax.axis("off")
        ax.set_title(title, fontsize=11.5, fontweight="bold", color=colour, loc="left", pad=12)
        for k, (cx, cy) in enumerate(centres):
            ax.add_patch(plt.Circle((cx, cy), radius, facecolor="#eef2f7",
                                    edgecolor=MUTED, lw=1.2, zorder=2))
            has_a, has_b = (True, True) if same else (k % 2 == 0, k % 2 == 1)
            # Place markers at fixed angles so the picture is reproducible and
            # the two panels are directly comparable.
            if has_a:
                for ang in (115, 155, 205, 245):
                    t = np.deg2rad(ang)
                    ax.plot(cx + radius * np.cos(t), cy + radius * np.sin(t), "o",
                            color=BLUE, ms=8, zorder=3)
            if has_b:
                for ang in (-25, 15, 55, 95):
                    t = np.deg2rad(ang)
                    ax.plot(cx + radius * np.cos(t), cy + radius * np.sin(t), "s",
                            color=ORANGE, ms=7.5, zorder=3)
        note = ("Every cell carries both, so the antibody can\nbridge them and recruit an immune cell."
                if same else
                "Averaged over the tumour this looks identical,\nbut no single cell carries both.")
        ax.text(5, -0.45, note, ha="center", va="top", fontsize=9.2,
                color=INK_2, linespacing=1.5)

    handles = [plt.Line2D([], [], marker="o", color=BLUE, ls="", ms=8, label="Target A"),
               plt.Line2D([], [], marker="s", color=ORANGE, ls="", ms=7.5, label="Target B")]
    leg = fig.legend(handles=handles, loc="upper right", frameon=False, fontsize=9.5,
                     bbox_to_anchor=(0.99, 1.02))
    for t in leg.get_texts():
        t.set_color(INK_2)
    fig.suptitle("Why we check 47,061 tumour cells one at a time", x=0.012, y=0.99,
                 ha="left", fontsize=13.5, fontweight="bold", color=INK)
    fig.tight_layout(rect=(0, 0.02, 1, 0.91))
    fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def fig_worked_example(gene: str = "TNFRSF17", label: str = "BCMA",
                       tumour_path=None, normal_path=None) -> dict:
    """One target, start to finish: where the tumour is, and what else carries it.

    BCMA is the ideal teaching case. It is in 100% of myeloma samples and
    essentially no other cancer, and its only substantial normal expression is
    in plasma cells - a lineage patients can live without, which is precisely
    why teclistamab is an approved drug despite that expression.
    """
    tumour_path = tumour_path or (C.FIGURES / "story_example_tumour.png")
    normal_path = normal_path or (C.FIGURES / "story_example_normal.png")

    _, frac = D.tumour_by_cancer_type()
    normal = D.load_normal_sc()
    pos = (frac.loc[gene] * 100).sort_values(ascending=False)
    nrm = normal.loc[gene].sort_values(ascending=False)

    # --- tumour side: emphasis on the one indication that matters ---
    top = pos.head(10).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    colours = [BLUE if v >= 50 else DEEMPH for v in top]
    ax.barh([c.replace("_", " ") for c in top.index], top.to_numpy(),
            color=colours, height=0.66)
    for y, v in enumerate(top.to_numpy()):
        ax.text(v + 1.6, y, f"{v:.0f}%", va="center", fontsize=9, color=INK_2)
    ax.set_xlim(0, 112)
    _style(ax, f"{label}: which cancers carry it",
           "Percent of that cancer's samples expressing the target",
           xlabel="% of tumour samples positive")
    ax.grid(axis="y", lw=0)
    fig.tight_layout(); fig.savefig(tumour_path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    # --- normal side: the single spike, and whether it is survivable ---
    topn = nrm.head(10).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.6, 4.4))
    vital = set(C.VITAL_CELL_TYPES)
    colours = [CRITICAL if c in vital else ORANGE if v > 100 else DEEMPH
               for c, v in zip(topn.index, topn.to_numpy())]
    ax.barh(list(topn.index), topn.to_numpy(), color=colours, height=0.66)
    for y, v in enumerate(topn.to_numpy()):
        ax.text(v + max(topn) * 0.02, y, f"{v:,.0f}", va="center", fontsize=9, color=INK_2)
    ax.set_xlim(0, max(topn) * 1.22)
    any_vital = any(c in vital for c in topn.index)
    sub = ("The safety question. Red marks a lineage whose loss would be dose-limiting"
           + ("." if any_vital else " - none appear here, which is the good news."))
    _style(ax, f"{label}: which healthy cells also carry it", sub,
           xlabel="expression in normal tissue (nCPM)")
    ax.grid(axis="y", lw=0)
    peak = nrm.index[0]
    ax.text(max(topn) * 0.42, 0.4,
            f"Only {peak} carry it.\nPatients can live without them -\nwhich is why this target became a drug.",
            fontsize=9, color=INK_2, linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.5", facecolor=SURFACE, edgecolor=GRID))
    fig.tight_layout(); fig.savefig(normal_path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)

    return {
        "gene": gene, "label": label,
        "top_cancer": pos.index[0].replace("_", " "),
        "top_cancer_pct": float(pos.iloc[0]),
        "second_cancer": pos.index[1].replace("_", " "),
        "second_cancer_pct": float(pos.iloc[1]),
        "normal_peak_cell": nrm.index[0],
        "normal_peak_value": float(nrm.iloc[0]),
        "normal_second_value": float(nrm.iloc[1]),
        "is_vital": nrm.index[0] in vital,
    }


def fig_scoreboard(validation: pd.DataFrame, cv: pd.DataFrame, path=None) -> None:
    """Simple rubric versus machine learning, on one axis.

    The headline of the project: the declared rubric beats the learned model,
    and the model falls below chance when asked about a pair type it was never
    shown. Chance is drawn, because every number here is only meaningful
    relative to it.
    """
    path = path or (C.FIGURES / "story_scoreboard.png")
    rows = []
    if validation is not None and not validation.empty:
        for grp, name in [("TAA|TAA [pathway_redundancy]", "Rubric\n(dual-target pairs)"),
                          ("TAA|checkpoint", "Rubric\n(tumour + checkpoint)"),
                          ("ALL", "Rubric\n(all pairs)")]:
            hit = validation[validation.group == grp]
            if not hit.empty and pd.notna(hit.iloc[0].get("roc_auc")):
                rows.append((name, float(hit.iloc[0]["roc_auc"]), BLUE))
    if cv is not None and not cv.empty:
        for key, name in [("arm-grouped CV", "Machine learning\n(fair test)"),
                          ("CD3E removed", "Machine learning\n(unseen target type)")]:
            hit = cv[cv.setting.astype(str).str.contains(key, regex=False)]
            if not hit.empty and pd.notna(hit.iloc[0].get("roc_auc")):
                rows.append((name, float(hit.iloc[0]["roc_auc"]), ORANGE))
    if not rows:
        return

    d = pd.DataFrame(rows, columns=["name", "auc", "colour"]).sort_values("auc")
    fig, ax = plt.subplots(figsize=(9.2, 0.72 * len(d) + 2.3))
    bars = ax.barh(d.name, d.auc, color=d.colour, height=0.62)
    ax.axvline(0.5, color=CRITICAL, lw=2, ls="--", zorder=3)
    ax.text(0.505, len(d) - 0.42, "coin flip (0.5)", color=CRITICAL, fontsize=9)
    for bar, v in zip(bars, d.auc):
        worse = v < 0.5
        ax.text(v + 0.012, bar.get_y() + bar.get_height() / 2,
                f"{v:.2f}" + ("  worse than guessing" if worse else ""),
                va="center", fontsize=9.5,
                color=CRITICAL if worse else INK_2,
                fontweight="bold" if worse else "normal")
    ax.set_xlim(0, 1.0)
    _style(ax, "The simple transparent method beat machine learning",
           "Higher is better. 0.5 means no better than a coin flip.",
           xlabel="accuracy at recognising real bispecific drugs (ROC-AUC)")
    ax.grid(axis="y", lw=0)
    handles = [mpatches.Patch(color=BLUE, label="Transparent rubric"),
               mpatches.Patch(color=ORANGE, label="Machine learning (LightGBM)")]
    # Below the axes: in-plot placement collides with the "worse than
    # guessing" annotation on the shortest bar.
    leg = ax.legend(handles=handles, frameon=False, fontsize=9.5, ncol=2,
                    loc="upper center", bbox_to_anchor=(0.5, -0.17))
    for t in leg.get_texts():
        t.set_color(INK_2)
    fig.tight_layout(); fig.savefig(path, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def build_all(universe_n: int, pairs_n: int, malignant_cells: int,
              validation: pd.DataFrame, cv: pd.DataFrame) -> dict:
    fig_pipeline(universe_n, pairs_n, malignant_cells)
    fig_same_cell_concept()
    fig_scoreboard(validation, cv)
    return fig_worked_example()
