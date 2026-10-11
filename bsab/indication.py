"""Per-indication target shortlists.

A single global ranking is the wrong output for this problem. A target's value
is indication-specific: CLDN18 matters in gastric cancer and nowhere else, BCMA
in myeloma, DLL3 in small-cell lung. The global rubric answers "is this a good
pair in general", which no development programme actually asks.

So for each cancer type we rank tumour antigens by positivity *within that
cancer type's cell lines*, penalised by normal-tissue exposure, and report the
shortlist. The safety terms are indication-independent (normal tissue does not
change), which is itself informative: the same antigen can be a good choice in
one indication and unusable in another purely because of how much tumour
carries it there.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C
from . import data_hpa as D


def rank_targets_by_indication(universe: pd.DataFrame, top_n: int = 15,
                               min_positivity: float = 0.30,
                               min_tumour_ntpm: float = 20.0) -> pd.DataFrame:
    """Shortlist tumour antigens per cancer type.

    `tumour_positivity` is the fraction of that cancer type's cell lines
    expressing the antigen at nTPM >= 10 - a prevalence estimate, i.e. roughly
    what fraction of patients might be addressable.

    **The ranking deliberately does not subtract normal-tissue expression.** An
    earlier version did, and it was wrong in an instructive way: subtracting
    measured normal expression rewards genes that are *poorly detected* just as
    much as genes that are genuinely tumour-restricted, so the shortlists filled
    with obscure proteins at near-zero normal signal while BCMA - 100% positive
    in myeloma lines - was pushed out for its (entirely expected) 566 nCPM in
    plasma cells. That expression in plasma cells is exactly the liability a
    triage needs to *see*, not something to silently penalise: plasma-cell
    depletion is the accepted trade for teclistamab.

    So antigens are ranked on tumour prevalence and magnitude among
    antibody-reachable proteins, and the normal-tissue burden is reported in
    adjacent columns for the reader to weigh. `normal_peak_cell_type` names the
    lineage carrying the most, which is usually the whole safety question.
    """
    taa = universe[universe.role == "TAA"].set_index("gene")
    mean_t, frac = D.tumour_by_cancer_type()
    frac = frac.reindex(taa.index).fillna(0.0)
    mean_t = mean_t.reindex(taa.index).fillna(0.0)
    normal = D.load_normal_sc().reindex(taa.index).fillna(0.0)

    vital = [c for c in C.VITAL_CELL_TYPES if c in normal.columns]
    normal_max = normal.max(axis=1)
    normal_peak_ct = normal.idxmax(axis=1)
    normal_vital = np.log1p(normal[vital].max(axis=1))
    breadth = (normal >= C.EXPR_EXPRESSED).sum(axis=1) / normal.shape[1]
    surface = (taa["surface_confidence"] if "surface_confidence" in taa.columns
               else pd.Series("unknown", index=taa.index))

    rows = []
    for ctype in frac.columns:
        pos = frac[ctype]
        level = mean_t[ctype]
        keep = (pos >= min_positivity) & (level >= min_tumour_ntpm) & surface.isin(
            ["high", "medium"])
        if not keep.any():
            continue
        g = pos.index[keep]
        sub = pd.DataFrame({
            "cancer_type": ctype,
            "gene": g,
            "tumour_positivity": pos[keep].to_numpy(),
            "tumour_mean_ntpm": level[keep].to_numpy(),
            "normal_peak_ncpm": normal_max[keep].to_numpy(),
            "normal_peak_cell_type": normal_peak_ct[keep].to_numpy(),
            "normal_vital_max_log": normal_vital[keep].to_numpy(),
            "normal_breadth": breadth[keep].to_numpy(),
            "surface_confidence": surface[keep].to_numpy(),
            "in_known_bsab": taa.loc[g, "in_known_bsab"].to_numpy(),
        })
        # Prevalence and magnitude, penalised by expression in *dose-limiting*
        # lineages and by how broadly the antigen is expressed at all - not by
        # its peak normal value. That distinction is the whole point: BCMA's
        # 566 nCPM in plasma cells is survivable (plasma-cell depletion is the
        # accepted trade for teclistamab), whereas B2M's 5,900 nCPM in
        # enterocytes is not. Peak-normal penalties cannot tell those apart and
        # rank housekeeping genes (B2M, GAPDH, ANXA2) at the top.
        sub["indication_score"] = (
            1.0 * sub.tumour_positivity
            + 0.25 * np.log1p(sub.tumour_mean_ntpm) / np.log1p(1000.0)
            + 0.10 * (sub.surface_confidence == "high").astype(float)
            - 0.16 * sub.normal_vital_max_log
            - 0.90 * sub.normal_breadth
        )
        sub = sub.sort_values(["indication_score", "tumour_mean_ntpm"],
                             ascending=False).head(top_n)
        sub.insert(0, "rank", range(1, len(sub) + 1))
        rows.append(sub)

    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    for c in ["tumour_positivity", "tumour_mean_ntpm", "normal_peak_ncpm",
              "normal_vital_max_log", "normal_breadth", "indication_score"]:
        out[c] = out[c].round(3)
    return out


def recovered_known_targets(by_ind: pd.DataFrame) -> pd.DataFrame:
    """Which documented-bispecific antigens each indication shortlist recovers.

    A weak but real sanity check: the shortlists should surface the targets
    already being developed in that indication.
    """
    if by_ind.empty:
        return pd.DataFrame()
    g = (by_ind.groupby("cancer_type")
         .agg(shortlisted=("gene", "size"),
              known_recovered=("in_known_bsab", "sum"),
              known_genes=("gene", lambda s: ", ".join(
                  by_ind.loc[s.index][by_ind.loc[s.index].in_known_bsab].gene)))
         .reset_index()
         .sort_values("known_recovered", ascending=False))
    return g
