"""Transparent, role-aware scoring of bispecific target pairs.

This is the project's primary ranking, deliberately *not* a learned model. With
64 documented pairs - 48 sharing CD3E as one arm - a supervised model mostly
learns "put CD3E on it", and its cross-validated numbers do not reflect
discovery performance. A declared rubric can be audited and argued with; see
`model.py` for the supervised comparison.

**Why the rubric is split by design class.** An earlier single rubric scored
*below chance* (ROC-AUC 0.43) because it rewarded same-cell co-expression for
every pair. That is the right criterion only for dual-tumour-antigen pairs. In a
T-cell engager the effector arm is supposed to be on T cells and *absent* from
tumour, so same-cell co-detection is near zero by design - the metric was
penalising the dominant class of real bispecifics. The arms are not
interchangeable, so the scoring cannot be symmetric in what it rewards:

* **TAA x effector** (T/NK-cell engager). Judged almost entirely on the tumour
  arm: how much tumour carries it, how little vital normal tissue does. The
  effector arm is a fixed, validated anchor and its own expression is intended.
* **TAA x TAA** (dual antigen). Here same-cell co-detection *is* the criterion,
  and the AND-gate minimum in normal tissue is what pairing buys.
* **TAA x checkpoint** (tumour-localised checkpoint blockade). Tumour-arm
  quality plus evidence the checkpoint is present in the microenvironment.
* **checkpoint x checkpoint / checkpoint x effector**. Expression data says very
  little about these; they are scored on a thin rubric and flagged low-confidence.

Because the classes are judged on different evidence, scores are percentiles
*within class* and the primary outputs are per-class rankings. A single global
ordering across classes would imply a comparability that does not exist.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Per class: component -> (weight, feature column, higher_is_better)
Rubric = dict[str, tuple[float, str, bool]]

TAA_EFFECTOR: Rubric = {
    "tumour_prevalence": (0.22, "taa_tumour_best_frac", True),
    "tumour_expression": (0.12, "taa_tumour_max_log", True),
    "breadth_across_cancers": (0.08, "taa_tumour_n_types_hit", True),
    "window_vs_vital_tissue": (0.18, "taa_window_vital_log", True),
    "vital_lineage_safety": (0.16, "taa_normal_vital_max_log", False),
    "normal_breadth_safety": (0.10, "taa_normal_frac_expressed", False),
    "lineage_specificity": (0.06, "taa_normal_tau", True),
    "surface_accessibility": (0.06, "taa_surface_confidence", True),
    "not_shed_to_blood": (0.02, "taa_secreted_to_blood", False),
}

# Dual-tumour-antigen bispecifics are built for two genuinely different
# reasons, so they get two rubrics and a pair surfaces if it is good under
# either. Blending them into one score buries both signals.
#
# (1) AND-gate selectivity / escape prevention: both antigens on the same
#     tumour cell, little normal tissue carrying both. CD19 x CD22 is the type.
TAA_TAA_ANDGATE: Rubric = {
    "same_cell_codetection": (0.24, "sc_pboth_max", True),
    "codetection_coincidence": (0.10, "sc_lift_max", True),
    "codetection_in_lines": (0.12, "line_double_pos_frac", True),
    "and_gate_normal_safety": (0.16, "normal_and_max_log", False),
    "vital_double_positive": (0.12, "normal_double_pos_vital_n", False),
    "and_gate_gain": (0.08, "and_gate_gain", True),
    "tme_selectivity": (0.08, "sc_tme_margin", True),
    "tumour_arm_prevalence": (0.04, "taa_tumour_best_frac", True),
    "surface_accessibility": (0.06, "taa_surface_confidence", True),
}

# (2) Pathway-redundancy blockade: hit two members of one signalling axis so the
#     tumour cannot route around either. EGFR x MET (amivantamab) and
#     HER2 x HER3 (zenocutuzumab) are the type. The decisive property is shared
#     function, not tumour selectivity - these targets are broadly expressed and
#     were never chosen for an expression window.
TAA_TAA_PATHWAY: Rubric = {
    "shared_molecular_function": (0.26, "mf_jaccard", True),
    "shared_biological_process": (0.16, "bp_jaccard", True),
    "coordinated_tumour_expression": (0.14, "tumour_profile_corr", True),
    "both_on_tumour": (0.12, "taa_tumour_best_frac", True),
    "tumour_expression_level": (0.10, "taa_tumour_max_log", True),
    "cancer_gene_involvement": (0.08, "cls_cancer_related__min", True),
    "druggable_precedent": (0.06, "taa_cls_fda_target", True),
    "vital_lineage_safety": (0.04, "taa_normal_vital_max_log", False),
    "surface_accessibility": (0.04, "taa_surface_confidence", True),
}

TAA_CHECKPOINT: Rubric = {
    "tumour_prevalence": (0.22, "taa_tumour_best_frac", True),
    "tumour_expression": (0.12, "taa_tumour_max_log", True),
    "window_vs_vital_tissue": (0.18, "taa_window_vital_log", True),
    "vital_lineage_safety": (0.16, "taa_normal_vital_max_log", False),
    "normal_breadth_safety": (0.10, "taa_normal_frac_expressed", False),
    "checkpoint_present_in_tme": (0.12, "sc_stromal_pboth_max", True),
    "surface_accessibility": (0.06, "taa_surface_confidence", True),
    "lineage_specificity": (0.04, "taa_normal_tau", True),
}

IMMUNE_ONLY: Rubric = {
    "expression_in_tme": (0.45, "sc_stromal_pboth_max", True),
    "codetection_coincidence": (0.25, "sc_lift_max", True),
    "normal_breadth_safety": (0.30, "normal_frac_expressed__max", False),
}

# class -> {rationale name: rubric}. A pair is scored under every rationale
# available to its class and keeps the best, with `rationale` recording which.
RUBRICS: dict[str, dict[str, Rubric]] = {
    "TAA|effector": {"effector_engager": TAA_EFFECTOR},
    "TAA|TAA": {"and_gate_selectivity": TAA_TAA_ANDGATE,
                "pathway_redundancy": TAA_TAA_PATHWAY},
    "TAA|checkpoint": {"tumour_localised_checkpoint": TAA_CHECKPOINT},
    "checkpoint|checkpoint": {"immune_combination": IMMUNE_ONLY},
    "checkpoint|effector": {"immune_combination": IMMUNE_ONLY},
}

# Classes where expression data is genuinely weak evidence.
LOW_CONFIDENCE_CLASSES = {"checkpoint|checkpoint", "checkpoint|effector"}


def _percentile(v: pd.Series, higher_is_better: bool) -> pd.Series:
    x = pd.to_numeric(v, errors="coerce")
    if not higher_is_better:
        x = -x
    return x.rank(pct=True, na_option="keep")


def score_pairs(features: pd.DataFrame) -> pd.DataFrame:
    """Score every pair under its design class's rubric.

    Percentiles are computed *within class*, so `score` ranks a pair against
    comparable designs rather than across incommensurable ones.
    """
    if "role_combo" not in features.columns:
        raise ValueError("features must carry 'role_combo'")

    out = pd.DataFrame(index=features.index)
    out["pair"] = features["pair"].to_numpy()
    out["design_class"] = features["role_combo"].to_numpy()
    out["score"] = np.nan
    out["coverage"] = np.nan
    out["rationale"] = None
    out["low_confidence"] = out.design_class.isin(LOW_CONFIDENCE_CLASSES)

    for cls, rationales in RUBRICS.items():
        rows = features.index[features["role_combo"] == cls]
        if len(rows) == 0:
            continue
        sub = features.loc[rows]
        for rname, rubric in rationales.items():
            total_w = sum(w for w, _, _ in rubric.values())
            wsum = pd.Series(0.0, index=rows)
            wpresent = pd.Series(0.0, index=rows)
            for cname, (w, col, hib) in rubric.items():
                colname = f"c_{rname}__{cname}"
                if colname not in out.columns:
                    out[colname] = np.nan
                if col not in sub.columns:
                    continue
                pct = _percentile(sub[col], hib)
                out.loc[rows, colname] = pct.to_numpy()
                avail = pct.notna()
                wsum[avail] += w * pct[avail]
                wpresent[avail] += w

            s = pd.Series(np.where(wpresent > 0, wsum / wpresent, np.nan), index=rows)
            cov = wpresent / total_w
            scol, ccol = f"score__{rname}", f"coverage__{rname}"
            for c in (scol, ccol):
                if c not in out.columns:
                    out[c] = np.nan
            out.loc[rows, scol] = s.to_numpy()
            out.loc[rows, ccol] = cov.to_numpy()

            # Keep the best rationale for this pair.
            better = s.notna() & (out.loc[rows, "score"].isna() | (s > out.loc[rows, "score"]))
            idx = rows[better.to_numpy()]
            out.loc[idx, "score"] = s[better].to_numpy()
            out.loc[idx, "coverage"] = cov[better].to_numpy()
            out.loc[idx, "rationale"] = rname
    return out


def validate_against_known(scored: pd.DataFrame, labels: pd.Series,
                           by_class: bool = True) -> pd.DataFrame:
    """Does the rubric rank documented pairs above chance, within each class?

    Reported as enrichment over the chance rate (k/n). With negatives
    unlabelled, enrichment and recall@k are the interpretable framings; ROC-AUC
    is included but is an upper bound, since many "negatives" are untried pairs.
    """
    from sklearn.metrics import roc_auc_score

    y = np.asarray(labels).astype(int)
    frames = []
    groups: list[tuple[str, np.ndarray, str]] = [("ALL", np.ones(len(scored), bool), "score")]
    if by_class:
        for cls in scored.design_class.unique():
            m = (scored.design_class == cls).to_numpy()
            groups.append((cls, m, "score"))
            # Where a class has several rationales, score each on its own so a
            # strong one is not hidden behind a weak one.
            for rname in RUBRICS.get(cls, {}):
                col = f"score__{rname}"
                if col in scored.columns and len(RUBRICS[cls]) > 1:
                    groups.append((f"{cls} [{rname}]", m, col))

    for name, mask, scol in groups:
        s = pd.to_numeric(scored.loc[mask, scol], errors="coerce").to_numpy()
        yy = y[mask]
        ok = ~np.isnan(s)
        s, yy = s[ok], yy[ok]
        n, npos = len(s), int(yy.sum())
        rec: dict[str, object] = {"group": name, "n_pairs": n, "n_known": npos}
        if npos and npos < n:
            rec["roc_auc"] = round(float(roc_auc_score(yy, s)), 4)
            ys = yy[np.argsort(-s)]
            for k in (50, 100, 500, 1000):
                if k > n:
                    continue
                hit = int(ys[:k].sum())
                rec[f"recall@{k}"] = round(hit / npos, 3)
                rec[f"enrich@{k}"] = round((hit / k) / (npos / n), 2)
            pcts = [float((s < v).mean()) for v in s[yy == 1]]
            rec["median_known_percentile"] = round(float(np.median(pcts)), 3)
        frames.append(rec)
    return pd.DataFrame(frames)


def component_table(rubrics: dict[str, dict[str, Rubric]] | None = None) -> pd.DataFrame:
    """The rubrics themselves, as a publishable table."""
    rubrics = rubrics or RUBRICS
    rows = []
    for cls, rationales in rubrics.items():
        for rname, rub in rationales.items():
            for name, (w, col, hib) in rub.items():
                rows.append({"design_class": cls, "rationale": rname, "component": name,
                             "weight": w, "feature": col, "higher_is_better": hib})
    return pd.DataFrame(rows)
