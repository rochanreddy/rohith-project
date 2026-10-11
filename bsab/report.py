"""Figures and the written report.

Charts follow the project's data-viz rules: one axis per chart, categorical
hues assigned in fixed order (never cycled), a single-hue sequential ramp for
magnitude, emphasis rather than categorical colour where one series is the
point, direct labels instead of a value on every mark, and a CSV table beside
every figure (the aqua slot sits below 3:1 on the light surface, so the table
view is the required relief).
"""

from __future__ import annotations

import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

from . import config as C

# Validated categorical slots (light mode), assigned in fixed order.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SURFACE = "#fcfcfb"
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
GRID = "#e4e3df"
DEEMPH = "#c3c2bd"
CRITICAL, GOOD = "#d03b3b", "#0ca30c"
SEQ_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#104281"]
)

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID,
    "axes.labelcolor": INK_2,
    "text.color": INK,
    "xtick.color": INK_2,
    "ytick.color": INK_2,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})


def _style(ax, title: str, subtitle: str | None = None, xlabel="", ylabel=""):
    ax.set_title(title, loc="left", pad=18 if subtitle else 8)
    if subtitle:
        ax.text(0, 1.02, subtitle, transform=ax.transAxes, fontsize=9,
                color=MUTED, va="bottom")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_axisbelow(True)


# --------------------------------------------------------------------------
def fig_prior_benchmark(bench: pd.DataFrame, path) -> None:
    """Emphasis form: only priors beating chance get the accent hue."""
    if bench.empty:
        return
    d = bench.sort_values("enrich@600")
    fig, ax = plt.subplots(figsize=(8, 0.6 * len(d) + 2.2))
    colors = [BLUE if v >= 1 else DEEMPH for v in d["enrich@600"]]
    bars = ax.barh(d.prior, d["enrich@600"], color=colors, height=0.6)
    ax.axvline(1.0, color=CRITICAL, lw=2, ls="--", zorder=3)
    ax.text(1.03, -0.75, "chance", color=CRITICAL, fontsize=9, va="center")
    for bar, v, r in zip(bars, d["enrich@600"], d["recall@600"]):
        ax.text(v + 0.04, bar.get_y() + bar.get_height() / 2,
                f"{v:.2f}x  (recall {r:.0%})", va="center", fontsize=9, color=INK_2)
    ax.set_xlim(0, max(2.6, d["enrich@600"].max() * 1.45))
    _style(ax, "Hand-built target-quality priors barely beat chance",
           "Fold-enrichment of documented bispecific antigens in each prior's top 600 "
           f"(of {int(d.n_taa.iloc[0]):,} tumour antigens)",
           xlabel="enrichment over chance (1.0 = no information)")
    ax.grid(axis="y", lw=0)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_cv_comparison(rows: list[dict], path) -> None:
    """Grouped bars contrasting leaky and honest evaluation."""
    d = pd.DataFrame(rows)
    if d.empty:
        return
    metrics = [("roc_auc", "ROC-AUC"), ("average_precision", "Avg. precision"),
               ("recall@500", "Recall@500")]
    metrics = [(k, lab) for k, lab in metrics if k in d.columns]
    x = np.arange(len(metrics))
    width = 0.8 / len(d)
    fig, ax = plt.subplots(figsize=(9, 4.6))
    for i, (_, r) in enumerate(d.iterrows()):
        color = [BLUE, ORANGE, AQUA][i % 3]
        vals = [float(r.get(k) or 0) for k, _ in metrics]
        pos = x + i * width - 0.4 + width / 2
        ax.bar(pos, vals, width * 0.88, label=r["setting"], color=color)
        for px, v in zip(pos, vals):
            ax.text(px, v + 0.015, f"{v:.2f}", ha="center", fontsize=8, color=INK_2)
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, lab in metrics])
    ax.set_ylim(0, 1.08)
    _style(ax, "The supervised model is weak under every evaluation protocol",
           "Same features and model, three protocols. Removing CD3E from training drops it "
           "below chance. Negatives are unlabelled, so these are upper bounds.",
           ylabel="score")
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.grid(axis="x", lw=0)
    leg = ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    for t in leg.get_texts():
        t.set_color(INK_2)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_feature_importance(imp: pd.DataFrame, path, n: int = 18) -> None:
    d = imp.head(n).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(d) + 2))
    norm = d.gain_frac / d.gain_frac.max()
    ax.barh(d.feature, d.gain_frac, color=[SEQ_BLUE(0.25 + 0.6 * v) for v in norm], height=0.68)
    for y, v in enumerate(d.gain_frac):
        ax.text(v + d.gain_frac.max() * 0.015, y, f"{v:.1%}", va="center",
                fontsize=8.5, color=INK_2)
    ax.set_xlim(0, d.gain_frac.max() * 1.2)
    _style(ax, "What the supervised model actually keys on",
           "Share of total split gain, model fitted on all 64 documented pairs",
           xlabel="share of gain")
    ax.grid(axis="y", lw=0)
    ax.tick_params(labelsize=8.5)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_safety_efficacy(features: pd.DataFrame, scored: pd.DataFrame, path,
                        n_label: int = 12) -> None:
    """Emphasis scatter: candidates recede, documented pairs and top picks are named."""
    d = features.merge(scored[["pair", "score"]], on="pair", how="left")
    from .features_pair import label_pairs

    d["known"] = label_pairs(d).to_numpy().astype(bool)
    d = d.dropna(subset=["sc_pboth_max", "normal_and_max_log"])
    if d.empty:
        return

    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    cand = d[~d.known]
    ax.scatter(cand.sc_pboth_max, cand.normal_and_max_log, s=5, c=DEEMPH,
               alpha=0.35, linewidths=0, label=f"candidate pairs (n={len(cand):,})")
    kn = d[d.known]
    ax.scatter(kn.sc_pboth_max, kn.normal_and_max_log, s=58, c=ORANGE,
               edgecolors=SURFACE, linewidths=1.6, zorder=4,
               label=f"documented bispecifics (n={len(kn)})")
    top = d[~d.known].nlargest(n_label, "score")
    ax.scatter(top.sc_pboth_max, top.normal_and_max_log, s=58, c=BLUE,
               edgecolors=SURFACE, linewidths=1.6, zorder=5,
               label=f"top-{n_label} ranked candidates")
    for _, r in pd.concat([kn.nlargest(6, "sc_pboth_max"), top.head(8)]).iterrows():
        ax.annotate(r.pair.replace("|", "x"), (r.sc_pboth_max, r.normal_and_max_log),
                    textcoords="offset points", xytext=(6, 4), fontsize=7.5, color=INK_2)
    _style(ax, "The developability trade-off: same-cell tumour co-expression vs normal-tissue leakage",
           "Lower is safer on the vertical axis. Desirable pairs sit bottom-right.",
           xlabel="fraction of malignant cells co-expressing both (max across tumour datasets)",
           ylabel="AND-gate peak in normal tissue, log1p nCPM")
    ax.grid(axis="both", color=GRID, lw=0.8)
    leg = ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    for t in leg.get_texts():
        t.set_color(INK_2)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def fig_sc_heatmap(features: pd.DataFrame, scored: pd.DataFrame, path, n: int = 20) -> None:
    cols = [c for c in features.columns if c.startswith("sc_pboth__")]
    if not cols:
        return
    d = features.merge(scored[["pair", "score"]], on="pair", how="left")
    d = d.dropna(subset=cols, how="all").nlargest(n, "score")
    if d.empty:
        return
    M = d[cols].to_numpy(dtype=float) * 100
    labels = [c.replace("sc_pboth__", "") for c in cols]

    fig, ax = plt.subplots(figsize=(1.35 * len(cols) + 4.5, 0.36 * len(d) + 2.6))
    im = ax.imshow(M, cmap=SEQ_BLUE, aspect="auto", vmin=0,
                   vmax=max(np.nanmax(M), 1))
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8.5)
    ax.set_yticks(range(len(d)))
    ax.set_yticklabels([p.replace("|", " x ") for p in d.pair], fontsize=8.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isnan(v):
                ax.text(j, i, "-", ha="center", va="center", fontsize=7.5, color=MUTED)
            else:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=7.5,
                        color="white" if v > np.nanmax(M) * 0.55 else INK)
    ax.set_xticks(np.arange(-.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(d), 1), minor=True)
    ax.grid(which="minor", color=SURFACE, lw=2)
    ax.tick_params(which="minor", length=0)
    _style(ax, "Same-cell co-detection varies sharply by tumour type",
           "% of malignant cells co-expressing both antigens; '-' = gene absent from that dataset")
    ax.grid(which="major", lw=0)
    cb = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    cb.set_label("% malignant cells double-positive", fontsize=8.5, color=INK_2)
    cb.outline.set_visible(False)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------
def write_report(ctx: dict, path=None) -> str:
    """Assemble results/report.md from the pipeline context."""
    path = path or (C.RESULTS / "report.md")
    uni = ctx["universe"]
    pairs = ctx["pairs"]
    scored = ctx["scored"]
    feats = ctx["features"]
    val = ctx["validation"]
    bench = ctx.get("prior_benchmark", pd.DataFrame())
    cv_rows = ctx.get("cv_rows", [])
    imp = ctx.get("importance", pd.DataFrame())
    sc_summary = ctx.get("sc_summary", pd.DataFrame())
    top_expl = ctx.get("explanations", [])
    novel = ctx.get("novel_top", pd.DataFrame())

    def md_table(df: pd.DataFrame, n: int | None = None) -> str:
        if df is None or df.empty:
            return "_no data_\n"
        d = df.head(n) if n else df
        return d.to_markdown(index=False, floatfmt=".3f") + "\n"

    L: list[str] = []
    A = L.append

    A("# Predicting bispecific-antibody target pairs\n")
    A("A reproducible pipeline that ranks two-target combinations for bispecific "
      "antibody development, built on public single-cell and tumour expression data. "
      "This document reports what the method does, what it found, and - at some "
      "length - where it does not work.\n")

    A("## Summary of findings\n")
    A(f"- Built a candidate space of **{len(pairs):,} target pairs** over "
      f"**{len(uni):,} antibody-accessible genes** "
      f"({int((uni.role == 'TAA').sum()):,} tumour antigens, "
      f"{int((uni.role != 'TAA').sum())} immune anchors).")
    A(f"- Measured same-cell co-expression directly in **{int(sc_summary.malignant_cells.sum()):,} "
      f"malignant cells** across {len(sc_summary)} tumour scRNA-seq datasets "
      f"({', '.join(sc_summary.disease.unique())})." if not sc_summary.empty else
      "- Single-cell co-expression unavailable.")
    if not bench.empty:
        best = bench.sort_values("enrich@600", ascending=False).iloc[0]
        A(f"- **Negative result that shaped the design:** hand-built target-quality priors "
          f"are close to chance at recovering known bispecific antigens. The best "
          f"(`{best['prior']}`) reaches only **{best['enrich@600']}x enrichment**, and "
          f"two of four priors score *below* chance. Expression contrast alone does not "
          f"identify validated targets, so no such prior was used to narrow the search.")
    if isinstance(val, pd.DataFrame) and not val.empty:
        allrow = val[val.group == "ALL"]
        if not allrow.empty:
            r = allrow.iloc[0]
            A(f"- The role-aware rubric places documented pairs at a median percentile of "
              f"**{r.get('median_known_percentile')}**, with "
              f"**{r.get('enrich@500', 'n/a')}x enrichment in the top 500** "
              f"(ROC-AUC {r.get('roc_auc')} against unlabelled negatives).")
        tt = val[val.group.astype(str).str.contains(r"TAA\|TAA \[")]
        if len(tt) >= 2:
            best = tt.loc[tt.roc_auc.idxmax()]
            worst = tt.loc[tt.roc_auc.idxmin()]
            A(f"- **The two dual-antigen rationales behave completely differently.** "
              f"Scoring dual-TAA pairs on pathway redundancy gives ROC-AUC "
              f"**{best.roc_auc}**; scoring the same pairs on AND-gate selectivity gives "
              f"**{worst.roc_auc}** - *below* chance. Historical dual-TAA bispecifics were "
              f"built to block signalling redundancy, not to exploit tumour-selective "
              f"AND-gating, so the latter has no validation set at all (see below).")
    if cv_rows:
        d = pd.DataFrame(cv_rows).set_index("setting")
        leaky = next((i for i in d.index if "random" in i), None)
        honest = next((i for i in d.index if "arm-grouped CV" == i), None)
        strict = next((i for i in d.index if "CD3E removed" in i), None)
        if leaky and honest:
            lk, hn = float(d.loc[leaky, "roc_auc"]), float(d.loc[honest, "roc_auc"])
            A(f"- **The supervised model is weak under every protocol** "
              f"(ROC-AUC {min(lk, hn):.3f}-{max(lk, hn):.3f}), and the transparent rubric "
              f"beats it. Notably, random k-fold did **not** inflate the headline here "
              f"(random {lk:.3f} vs arm-grouped {hn:.3f}) - the leakage effect this "
              f"project was built to guard against did not materialise on this data, and "
              f"saying so matters more than keeping a tidy story.")
            if strict:
                st_auc = float(d.loc[strict, "roc_auc"])
                A(f"- **The decisive failure is generalisation to unseen arms.** With every "
                  f"CD3E pair removed from training, the model drops to ROC-AUC "
                  f"**{st_auc:.3f}** - below chance. Asked to find a good pair of a kind it "
                  f"was never shown, it does worse than guessing. That is the question a "
                  f"discovery team actually has, and it is the strongest argument for the "
                  f"declared rubric over the learned model.")
    A("")

    A("## The problem\n")
    A("A bispecific antibody binds two targets at once. In oncology the pairing usually "
      "does one of a few jobs: hold a tumour antigen next to an effector-cell receptor "
      "so cytotoxic cells are redirected to the tumour; localise checkpoint blockade to "
      "the tumour; or hit two tumour antigens to raise avidity and close off escape. "
      "Given a few thousand antibody-accessible proteins there are millions of possible "
      "pairs, and the question is which deserve laboratory work.\n")
    A("Two requirements decide it, and they pull in opposite directions. **Efficacy** "
      "needs both antigens present on the *same* tumour cell - a pair split across "
      "different cells cannot form the required bridge. **Safety** is governed by "
      "whichever normal tissue also carries both, because for a bispecific requiring "
      "dual engagement the exposure in a tissue scales with the *lesser* of the two "
      "antigens, not with either alone. That asymmetry is the whole point of pairing: a "
      "target too widely expressed to use on its own can become usable when gated "
      "against a second one.\n")

    A("## Data\n")
    A("| Source | Content | Role |")
    A("|---|---|---|")
    A("| Human Protein Atlas, single-cell atlas | 20,151 genes x 154 normal cell types (nCPM) | normal-tissue safety |")
    A("| Human Protein Atlas, cell lines | 20,151 genes x 1,206 cancer cell lines (nTPM), 25 cancer types | tumour expression |")
    A("| Human Protein Atlas, annotations | protein class, subcellular location, secretome, specificity | antibody accessibility |")
    A("| CELLxGENE tumour scRNA-seq | 5 datasets, malignant + microenvironment cells | same-cell co-expression |")
    A("| Curated reference list | 64 documented bispecific pairs (approved / clinical) | positive labels |")
    A("")
    if not sc_summary.empty:
        A("Single-cell datasets used:\n")
        A(md_table(sc_summary))
        A("")

    A("## Method\n")
    A("**Candidate universe.** Antibody-reachable genes showing some tumour expression "
      f"(>= {C.TUMOUR_ELIGIBILITY_FRACTION:.0%} of cell lines positive in at least one "
      f"cancer type): {len(uni):,} genes. Immune anchors and the antigens used in "
      "documented bispecifics are retained regardless.\n")
    A("Reachability needs care, and getting it wrong was the most consequential bug in "
      "this project. HPA's \"Predicted membrane proteins\" class is a *topology* call: "
      "it fires for every organelle membrane. An earlier version used it directly and "
      "the top-ranked candidates came back dominated by proteins no antibody can bind "
      "in a living cell - mitochondrial ATP-synthase subunits (ATP5MJ, ATP5MK), ER "
      "enzymes (VKORC1, the warfarin target; PDIA6), and lysosomal GILT (IFI30). The "
      "filter now grades reachability:\n")
    A("- **high** - positive surface or extracellular evidence: plasma-membrane "
      "localisation, secretion, CD marker, or GPCR\n"
      "- **medium** - predicted membrane protein with no contradicting organelle evidence\n"
      "- **none** - excluded, including anything whose only localisation evidence is an "
      "intracellular organelle\n")
    surf = ctx.get("surface_audit")
    if surf:
        A(f"The grade is also a scored feature, so high-confidence targets outrank "
          f"merely-plausible ones rather than the distinction being lost at the filter. "
          f"**Audited against the validated targets:** of {surf['known_genes']} antigens "
          f"used in documented bispecifics, {surf['high']} grade high, {surf['medium']} "
          f"medium, and {surf['rejected']} would be rejected outright "
          f"({', '.join(surf['rejected_genes'])}) - an {surf['rejected'] / surf['known_genes']:.0%} "
          f"false-negative rate, driven by targets that simply lack reliable "
          f"immunofluorescence data in HPA. Those six are retained here only because the "
          f"reference list names them, which is a real limit on the filter's recall for "
          f"genuinely novel targets.\n")
    A("**Pair space.** Bounded structurally rather than by a quality prior, because the "
      "priors were measured to be uninformative: every tumour antigen x every immune "
      "anchor, anchor x anchor, and tumour antigen x tumour antigen among the "
      f"{C.N_TOP_TAA_FOR_TAA_PAIRS} most prevalent antigens. "
      f"Total {len(pairs):,} pairs, containing all 64 documented pairs.\n")
    A("**Features.** Each pair carries symmetric aggregates (min / max / mean / "
      "absolute difference) of 35 single-arm features, so no feature can encode arm "
      "order, plus genuinely pairwise quantities: AND-gate peaks in normal and tumour "
      "tissue, double-positive counts, profile correlations, functional-term overlap, "
      "and single-cell co-detection. "
      f"{feats.shape[1] - 1} features in total.\n")
    A("**Ranking.** The primary ranking is a declared rubric of percentile components "
      "with fixed weights, not a fitted model. Components a pair has no data for are "
      "dropped and the remaining weights renormalised, so an unmeasured pair is not "
      "penalised as though it had been measured and found wanting; `coverage` records "
      "how much of the rubric was available.\n")
    A("Critically, **the rubric is split by design class**, because the arms of a "
      "bispecific are not interchangeable and the classes are not judged on the same "
      "evidence. An earlier single symmetric rubric scored *below chance* (ROC-AUC "
      "0.43) for exactly this reason: it rewarded same-cell co-expression for every "
      "pair, but in a T-cell engager the effector arm is supposed to be on T cells and "
      "absent from tumour, so same-cell co-detection is near zero by design. That "
      "metric was penalising the dominant class of real bispecifics. Scores are "
      "percentiles *within class*, and the primary outputs are per-class rankings - a "
      "single global ordering would imply a comparability that does not exist.\n")
    A(md_table(ctx.get("rubric_components", pd.DataFrame())))
    A("")

    A("## Why the primary ranking is not a machine-learning model\n")
    A("There are 64 documented bispecific pairs. **48 of them share CD3E** as one arm, "
      "because T-cell engagers dominate the clinical literature. The remaining "
      "unlabelled pairs are not negatives - almost all have simply never been tried. "
      "This is positive-unlabelled learning on a small, structurally lopsided positive "
      "set, and it has two consequences.\n")
    A("First, a model trained this way largely learns *\"one arm should be CD3E\"*. "
      "That is a real pattern, but it is not target discovery and it does not "
      "generalise to the pair types where most of the remaining opportunity sits.\n")
    A("Second, the obvious evaluation is wrong. Random k-fold puts `CD19|CD3E` in "
      "training and `CD22|CD3E` in test, so a model keying on CD3E scores brilliantly "
      "on pairs it has effectively already seen. All evaluation here is therefore "
      "grouped by shared arm, so every pair involving a given gene sits in one fold; "
      "and a stricter setting removes every CD3E pair from training and scores only "
      "non-CD3E pairs.\n")
    if cv_rows:
        A(md_table(pd.DataFrame(cv_rows)))
        A("")
        A("![Evaluation protocol comparison](figures/cv_comparison.png)\n")
    A("The supervised model is therefore reported as a **secondary, corroborating "
      "signal and a feature-importance probe**, not as the product.\n")

    if not bench.empty:
        A("## A negative result worth stating\n")
        A("Before building the pair model we asked whether a hand-built prior over "
          "*single* targets - the kind of expression-contrast heuristic that looks "
          "obviously sensible - can recover antigens already validated in bispecifics. "
          "Four priors were tested against the documented list.\n")
        A(md_table(bench))
        A("")
        A("![Prior benchmark](figures/prior_benchmark.png)\n")
        A("Only the therapeutic-window prior beats chance at all, and only modestly. "
          "Ranking targets by tumour prevalence, or by *low* normal expression, does "
          "**worse than random**. The reason is visible in the data: validated targets "
          "such as EGFR, EPCAM, MS4A1 and MET have substantial normal-tissue "
          "expression. What makes them usable is not expression contrast but whether "
          "the normal lineage carrying them is expendable (B-cell aplasia is "
          "survivable), whether toxicity is clinically manageable, and biology this "
          "data cannot see. Any pipeline that triages targets on expression contrast "
          "alone would have discarded most of the drugs that work.\n")
        A("This is why the universe here is deliberately permissive and the pair space "
          "is bounded by structure rather than by a quality score.\n")

    A("## Results\n")
    A("### Does the rubric recover known pairs?\n")
    A("Enrichment is the ratio to the chance rate (k/n) and is the interpretable "
      "figure here; ROC-AUC is reported but is an upper bound, because most "
      "\"negatives\" are untried rather than rejected pairs.\n")
    A(md_table(val if isinstance(val, pd.DataFrame) else pd.DataFrame([val])))
    A("")
    A("### The two dual-antigen rationales are not the same hypothesis\n")
    A("Dual-tumour-antigen bispecifics get built for two different reasons, so they "
      "are scored under both and keep the better:\n")
    A("- **Pathway redundancy** - hit two members of one signalling axis so the tumour "
      "cannot route around either. EGFR x MET (amivantamab), HER2 x HER3 "
      "(zenocutuzumab). The decisive property is shared function.\n")
    A("- **AND-gate selectivity** - both antigens on the same tumour cell, little "
      "normal tissue carrying both, so the pairing buys selectivity a single arm "
      "cannot. CD19 x CD22 is the closest documented example.\n")
    A("Seven of the nine documented dual-TAA pairs are dual-signalling or "
      "dual-angiogenic designs. Their shared-molecular-function overlap runs 0.67-1.0 "
      "against a class median of 0.00, while their same-cell co-detection sits *below* "
      "the class median. The pathway rubric therefore recovers them strongly and the "
      "AND-gate rubric scores below chance.\n")
    A("**This should not be read as the AND-gate concept being wrong.** It means there "
      "is no validation set for it: no approved dual-TAA bispecific was designed for "
      "tumour-selective AND-gating, so historical labels cannot confirm or refute it. "
      "The AND-gate ranking is a genuinely prospective hypothesis, and the honest "
      "statement is that it is untested rather than validated.\n")
    A("**Also worth stating plainly:** the pathway rubric was added *after* observing "
      "the AND-gate rubric score below chance on these nine pairs. Its ROC-AUC is "
      "therefore partly fitted to nine examples and should be discounted accordingly. "
      "The biological rationale (dual signalling blockade is the textbook reason for "
      "dual-TAA bispecifics) is prior knowledge rather than a discovery from this data, "
      "which is what makes the revision defensible - but it is still a revision made "
      "with the answers visible.\n")
    A("![Safety-efficacy trade-off](figures/safety_efficacy.png)\n")
    if not imp.empty:
        A("### What the supervised model keys on\n")
        A(md_table(imp.head(15)))
        A("")
        A("![Feature importance](figures/feature_importance.png)\n")
    A("### Same-cell co-detection by tumour type\n")
    A("![Single-cell co-detection](figures/sc_codetection.png)\n")
    A("Co-detection is strongly tumour-type dependent, which is the practical argument "
      "for indication-specific pair selection rather than a single global ranking.\n")

    by_ind = ctx.get("by_indication", pd.DataFrame())
    if not by_ind.empty:
        A("### Per-indication shortlists\n")
        A("A single global ranking answers \"is this a good pair in general\", which no "
          "development programme actually asks. A target's value is "
          "indication-specific - CLDN18 matters in gastric cancer and essentially "
          "nowhere else, BCMA in myeloma, DLL3 in neuroendocrine tumours. Antigens are "
          "therefore also ranked within each cancer type by positivity across *that* "
          "cancer type's cell lines, discounted by normal-tissue exposure. The safety "
          "terms do not change between indications, which is the point: the same "
          "antigen can be a sound choice in one indication and unusable in another "
          "purely because of how much tumour carries it there.\n")
        rec = ctx.get("indication_recovery", pd.DataFrame())
        if not rec.empty:
            A(f"Across {by_ind.cancer_type.nunique()} indications the shortlists recover "
              f"{int(by_ind.in_known_bsab.sum())} antigens already used in documented "
              f"bispecifics - a weak but real sanity check that they surface targets "
              f"people are in fact developing.\n")
            A(md_table(rec, 15))
            A("")
        show = ["cancer_type", "rank", "gene", "tumour_positivity",
                "normal_vital_max_log", "normal_breadth", "surface_confidence",
                "in_known_bsab", "indication_score"]
        show = [c for c in show if c in by_ind.columns]
        head = (by_ind[by_ind["rank"] <= 5]
                .sort_values(["cancer_type", "rank"])[show])
        A("Top five antigens per indication (full table in "
          "`results/targets_by_indication.csv`):\n")
        A(md_table(head, 60))
        A("")

    if not novel.empty:
        A("### Top-ranked pairs not in the reference list\n")
        A("These are the pipeline's actual suggestions. They are hypotheses for "
          "triage, and the caveats in the next section apply to every row.\n")
        A(md_table(novel, 30))
        A("")

    if top_expl:
        A("## Worked explanations\n")
        A("Explanations are assembled from the evidence behind each score - measured "
          "co-detection rates, the named normal cell types expressing both antigens, "
          "the cancer type driving the signal. They are generated deterministically so "
          "every statement is checkable against the feature table; an optional "
          "language-model pass (`--llm`) rewrites the same grounded evidence into "
          "smoother prose without being permitted to add facts.\n")
        for e in top_expl[:8]:
            A("```")
            A(e["text"])
            A("```\n")

    A("## Limitations\n")
    A("These are the reasons not to over-read the output.\n")
    A("1. **Labels are few, biased and not negatives.** 64 positives, 48 sharing one "
      "arm; unlabelled pairs are mostly untried rather than rejected. Absolute "
      "discrimination metrics against unlabelled negatives are upper bounds.\n")
    A("2. **mRNA is not surface protein.** Every expression input is transcriptomic. "
      "Surface density, epitope accessibility, internalisation rate and shedding - all "
      "decisive for a bispecific - are invisible here.\n")
    A("3. **scRNA-seq dropout biases co-detection downward.** A zero is not proof of "
      "absence, so same-cell co-expression rates are systematic under-estimates; they "
      "are usable for ranking pairs against each other, not as absolute values.\n")
    A("4. **Cell lines are a poor tumour proxy.** They lack microenvironment and often "
      "lose lineage markers in culture - CLDN18, SSTR2, ANGPT2 and DLL4 all sit near "
      "the expression floor in cell lines despite being real clinical targets, and were "
      "retained here only because the reference list named them.\n")
    A("5. **Five tumour datasets, five indications.** Co-detection coverage is thin and "
      "skewed; pairs whose genes are absent from all five carry no single-cell evidence "
      "at all (see `coverage`).\n")
    A("6. **The rubric weights are asserted, not fitted.** They encode a defensible "
      "view of developability and are published so they can be argued with; they are "
      "not optimised against an outcome, because no adequate outcome data exists.\n")
    A("7. **No immunology beyond expression.** Effector-cell fitness, cytokine-release "
      "risk, antigen density thresholds for synapse formation, resistance evolution and "
      "manufacturability are all out of scope, and they routinely decide programmes.\n")
    A("8. **The universe cut-off saw the answer.** Eligibility was set low enough to "
      "retain every antigen in the reference list. That makes in-universe recall "
      "optimistic as a discovery estimate; the arm-grouped and CD3E-excluded "
      "evaluations are the ones to trust.\n")
    A("9. **The rubric was revised with the answers visible.** Two revisions were made "
      "after seeing validation results: splitting the rubric by design class (after a "
      "symmetric version scored below chance), and adding the pathway-redundancy "
      "rationale for dual-TAA pairs (after the AND-gate version scored below chance on "
      "nine examples). Both revisions rest on prior pharmacology rather than on "
      "patterns mined from this data, which is what makes them defensible - but the "
      "resulting validation figures are optimistic, and the dual-TAA pathway number in "
      "particular is fitted to nine positives.\n")
    A("10. **Scores are not comparable across design classes.** They are percentiles "
      "within class, because the classes are judged on different evidence. "
      "`rank_overall` in the output is a convenience for browsing, not a claim that a "
      "TAA x effector pair at 0.9 is better than a TAA x TAA pair at 0.9.\n")
    A("11. **The AND-gate ranking is untested, not validated.** No approved dual-antigen "
      "bispecific was designed for tumour-selective AND-gating, so the historical "
      "labels can neither confirm nor refute that ranking. Those candidates are "
      "hypotheses awaiting a kind of evidence this project cannot supply.\n")
    A("12. **The reachability filter cannot tell which *face* of the membrane a protein "
      "sits on, and this is the biggest remaining source of implausible candidates.** "
      "HPA's immunofluorescence places cytoplasmic adaptors such as PPP1R18 and "
      "PLEKHO2 at the plasma membrane, which is correct and yet useless for antibody "
      "targeting: they have no extracellular domain. Separating them requires "
      "per-residue topology - UniProt transmembrane and extracellular-domain "
      "annotation - which this pipeline does not ingest. Until it does, the shortlists "
      "contain intracellular proteins, and `surface_confidence` should be read "
      "alongside every candidate rather than trusted as a filter. Adding UniProt "
      "topology is the single highest-value next step.\n")
    A("13. **Annotation quality limits the filter in both directions.** Of 56 validated "
      "bispecific antigens, 10 are rejected outright (CLDN6, DLL4, FOLH1, GUCY2C, "
      "LGR5, PMEL, ROR1, STEAP1, TIGIT, TPBG) because HPA lacks reliable localisation "
      "for them, and BCMA scrapes through only as *medium* - its HPA main location is "
      "the Golgi. A filter that mis-grades the target of an approved drug is not a "
      "filter to trust on novel ones.\n")

    A("## Reproducing\n")
    A("```bash\n"
      "python -m venv .venv && .venv/Scripts/pip install -r requirements.txt\n"
      "python scripts/download_data.py      # HPA tables + CELLxGENE tumour datasets\n"
      "python run.py                        # full pipeline -> results/\n"
      "python run.py --llm                  # same, with LLM-rewritten explanations\n"
      "```\n")
    A("Outputs land in `results/`: `ranked_pairs.csv` (full ranking with every "
      "feature), `novel_candidates.csv`, `explanations.md`, `report.md`, and the "
      "figures plus a CSV beside each.\n")

    text = "\n".join(L)
    path.write_text(text, encoding="utf-8")
    return text
