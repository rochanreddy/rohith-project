"""Per-target (single-arm) features.

These describe one candidate target in isolation: is it reachable by an
antibody, how much tumour carries it, and how much *normal* tissue carries it.
The last point is what kills most oncology targets, so normal-tissue exposure
is measured several ways (peak, breadth, and peak restricted to the lineages
whose loss is dose-limiting).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C
from . import data_hpa as D

# HPA protein-class flags worth carrying as features.
CLASS_FLAGS = {
    "cls_membrane": "Predicted membrane proteins",
    "cls_secreted": "Predicted secreted proteins",
    "cls_intracellular": "Predicted intracellular proteins",
    "cls_cd_marker": "CD markers",
    "cls_fda_target": "FDA approved drug targets",
    "cls_potential_target": "Potential drug targets",
    "cls_cancer_related": "Cancer-related genes",
    "cls_disease_related": "Disease related genes",
    "cls_gpcr": "G-protein coupled receptors",
    "cls_transporter": "Transporters",
    "cls_enzyme": "Enzymes",
    "cls_plasma_protein": "Plasma proteins",
    "cls_transcription_factor": "Transcription factors",
    "cls_ribosomal": "Ribosomal proteins",
}


def _log1p(s: pd.Series) -> pd.Series:
    return np.log1p(s.clip(lower=0))


def build(universe: pd.DataFrame) -> pd.DataFrame:
    genes = list(universe.gene)
    annot = D.load_annotations().set_index("gene").reindex(genes)
    normal = D.load_normal_sc().reindex(genes)
    mean_t, frac_t = D.tumour_by_cancer_type()
    mean_t = mean_t.reindex(genes)
    frac_t = frac_t.reindex(genes)
    subcell = D.load_subcellular().set_index("gene").reindex(genes)

    f = pd.DataFrame(index=pd.Index(genes, name="gene"))

    # ---- annotation / druggability ----
    pclass = annot["Protein class"].map(D.split_terms)
    for col, label in CLASS_FLAGS.items():
        f[col] = pclass.map(lambda s, lab=label: lab in s).astype(float)

    main_loc = subcell["Main location"].fillna("")
    f["loc_plasma_membrane"] = main_loc.str.contains("Plasma membrane", case=False).astype(float)
    f["loc_nuclear"] = main_loc.str.contains("Nucle", case=False).astype(float)

    # Graded antibody reachability from targets._accessibility: 1.0 = positive
    # surface/extracellular evidence, 0.5 = predicted membrane with nothing
    # contradicting it. Carried as a feature so high-confidence targets outrank
    # merely-plausible ones instead of the distinction being lost at the filter.
    uni_idx = universe.set_index("gene")
    f["surface_confidence"] = (
        uni_idx["surface_confidence"].reindex(genes)
        .map({"high": 1.0, "medium": 0.5}).fillna(0.0)
    )
    f["secreted_to_blood"] = (
        annot["Secretome location"].fillna("").str.contains("blood", case=False).astype(float)
    )
    f["n_interactions_log"] = _log1p(pd.to_numeric(annot["Interactions"], errors="coerce").fillna(0))

    for src, dst in [
        ("RNA tissue specificity score", "spec_tissue"),
        ("RNA single cell type specificity score", "spec_singlecell"),
        ("RNA cancer specificity score", "spec_cancer"),
    ]:
        f[dst] = pd.to_numeric(annot[src], errors="coerce")

    # ---- normal-tissue exposure (safety) ----
    nz = normal.fillna(0.0)
    f["normal_max_log"] = _log1p(nz.max(axis=1))
    f["normal_mean_log"] = _log1p(nz.mean(axis=1))
    f["normal_n_expressed"] = (nz >= C.EXPR_EXPRESSED).sum(axis=1).astype(float)
    f["normal_frac_expressed"] = f["normal_n_expressed"] / nz.shape[1]
    f["normal_tau"] = nz.apply(lambda r: D.tau_specificity(r.to_numpy()), axis=1)

    vital = [c for c in C.VITAL_CELL_TYPES if c in nz.columns]
    f["normal_vital_max_log"] = _log1p(nz[vital].max(axis=1))
    f["normal_vital_n_expressed"] = (nz[vital] >= C.EXPR_EXPRESSED).sum(axis=1).astype(float)

    # ---- tumour expression (efficacy) ----
    mt = mean_t.fillna(0.0)
    ft = frac_t.fillna(0.0)
    f["tumour_max_log"] = _log1p(mt.max(axis=1))
    f["tumour_mean_log"] = _log1p(mt.mean(axis=1))
    f["tumour_best_frac"] = ft.max(axis=1)
    f["tumour_mean_frac"] = ft.mean(axis=1)
    f["tumour_n_types_hit"] = (ft >= C.MIN_TUMOUR_LINE_FRACTION).sum(axis=1).astype(float)

    # ---- therapeutic window ----
    # Positive = tumour expression exceeds the worst normal exposure.
    f["window_log"] = f["tumour_max_log"] - f["normal_max_log"]
    f["window_vital_log"] = f["tumour_max_log"] - f["normal_vital_max_log"]

    f["role"] = universe.set_index("gene")["role"].reindex(genes)
    return f


def normal_profiles(universe: pd.DataFrame) -> pd.DataFrame:
    """log1p normal-cell-type expression matrix for the universe (pair features)."""
    normal = D.load_normal_sc().reindex(list(universe.gene)).fillna(0.0)
    return np.log1p(normal.clip(lower=0))


def tumour_profiles(universe: pd.DataFrame) -> pd.DataFrame:
    """log1p per-cancer-type mean tumour expression for the universe."""
    mean_t, _ = D.tumour_by_cancer_type()
    mt = mean_t.reindex(list(universe.gene)).fillna(0.0)
    return np.log1p(mt.clip(lower=0))
