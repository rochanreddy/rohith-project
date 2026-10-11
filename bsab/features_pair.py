"""Pair-level feature construction.

Two kinds of feature live here:

1. Symmetric aggregates of the two single-arm feature vectors. A bispecific has
   no intrinsic arm order, so every feature must be invariant to swapping the
   arms: min / max / mean / |difference| rather than (a, b) side by side. Raw
   per-arm columns would let a model memorise arm order and inflate
   cross-validation scores.

2. Genuinely pairwise quantities that cannot be derived from either arm alone:
   co-expression in tumour, co-expression in normal tissue, and the selectivity
   an AND-gate actually buys.

The central safety idea is the AND-gate. For a bispecific that needs both
antigens on the same cell, exposure in a normal cell type scales with
min(expr_a, expr_b), not with either arm alone. `normal_and_max_log` is the peak
of that minimum over all 154 normal cell types, and `and_gate_gain` is how far
it sits below the worse single arm - i.e. the safety the pairing buys.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C
from . import data_hpa as D
from . import pairops as P

AGG_SKIP = {"role"}

# Single-arm features worth resolving *to the tumour arm specifically*, with
# whether a high value is desirable. For a TAA x effector bispecific the two
# arms are not interchangeable: the effector arm (CD3E) is supposed to be on
# T cells and absent from tumour, so aggregates like min/max over both arms
# blur the quantity that actually decides the design - the tumour arm's own
# prevalence and normal-tissue exposure.
TAA_ARM_SPEC: dict[str, bool] = {
    "tumour_best_frac": True,
    "tumour_max_log": True,
    "tumour_n_types_hit": True,
    "window_log": True,
    "window_vital_log": True,
    "normal_max_log": False,
    "normal_vital_max_log": False,
    "normal_frac_expressed": False,
    "normal_vital_n_expressed": False,
    "normal_tau": True,
    "spec_cancer": True,
    "spec_singlecell": True,
    "loc_plasma_membrane": True,
    "surface_confidence": True,
    "secreted_to_blood": False,
    "cls_cd_marker": True,
    "cls_fda_target": True,
}


def _taa_arm_features(single: pd.DataFrame, ia: np.ndarray, ib: np.ndarray,
                      role: np.ndarray) -> pd.DataFrame:
    """Resolve selected features to the tumour-antigen arm.

    One TAA arm -> that arm's value. Two TAA arms -> the conservative one (the
    worse of the two), so a dual-antigen pair is judged on its weaker half
    rather than flattered by its stronger. No TAA arm -> NaN.
    """
    taa_a = role[ia] == "TAA"
    taa_b = role[ib] == "TAA"
    both = taa_a & taa_b
    only_a = taa_a & ~taa_b
    only_b = taa_b & ~taa_a

    out = {"n_taa_arms": (taa_a.astype(np.int8) + taa_b.astype(np.int8)).astype(np.float32)}
    for col, higher_is_better in TAA_ARM_SPEC.items():
        if col not in single.columns:
            continue
        v = single[col].to_numpy(dtype=np.float32)
        va, vb = v[ia], v[ib]
        conservative = np.minimum(va, vb) if higher_is_better else np.maximum(va, vb)
        vals = np.full(len(ia), np.nan, dtype=np.float32)
        vals[both] = conservative[both]
        vals[only_a] = va[only_a]
        vals[only_b] = vb[only_b]
        out[f"taa_{col}"] = vals
    return pd.DataFrame(out)


def _symmetric_aggregates(single: pd.DataFrame, ia: np.ndarray, ib: np.ndarray) -> pd.DataFrame:
    cols = [c for c in single.columns if c not in AGG_SKIP]
    M = single[cols].astype("float32").to_numpy()
    out = {}
    for j, c in enumerate(cols):
        a, b = M[ia, j], M[ib, j]
        out[f"{c}__min"] = np.minimum(a, b)
        out[f"{c}__max"] = np.maximum(a, b)
        out[f"{c}__mean"] = (a + b) / 2.0
        out[f"{c}__absdiff"] = np.abs(a - b)
    return pd.DataFrame(out)


def build(
    universe: pd.DataFrame,
    pairs: pd.DataFrame,
    single: pd.DataFrame,
    sc_features: pd.DataFrame | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    genes = list(universe.gene)
    gidx = {g: i for i, g in enumerate(genes)}
    ia = pairs.gene_a.map(gidx).to_numpy()
    ib = pairs.gene_b.map(gidx).to_numpy()

    if verbose:
        print(f"  symmetric aggregates over {len([c for c in single.columns if c not in AGG_SKIP])} single-arm features")
    feats = _symmetric_aggregates(single, ia, ib)
    feats.insert(0, "pair", pairs.pair.to_numpy())
    feats.insert(1, "gene_a", pairs.gene_a.to_numpy())
    feats.insert(2, "gene_b", pairs.gene_b.to_numpy())

    from .features_single import normal_profiles, tumour_profiles

    # ---------- normal tissue: AND-gate safety ----------
    if verbose:
        print("  normal-tissue co-expression (154 cell types)")
    normal_log = normal_profiles(universe)
    normal_cols = list(normal_log.columns)
    NL = normal_log.to_numpy(dtype=np.float32)
    NRAW = np.expm1(NL)

    feats["normal_and_max_log"] = np.log1p(P.max_of_min(NRAW, ia, ib))
    feats["normal_double_pos_n"] = P.both_above(NRAW, ia, ib, C.EXPR_EXPRESSED)
    feats["normal_either_pos_n"] = P.either_above(NRAW, ia, ib, C.EXPR_EXPRESSED)
    feats["normal_profile_corr"] = P.row_corr(NL, ia, ib)

    vital_idx = [i for i, c in enumerate(normal_cols) if c in C.VITAL_CELL_TYPES]
    if vital_idx:
        VR = np.ascontiguousarray(NRAW[:, vital_idx])
        feats["normal_and_vital_max_log"] = np.log1p(P.max_of_min(VR, ia, ib))
        feats["normal_double_pos_vital_n"] = P.both_above(VR, ia, ib, C.EXPR_EXPRESSED)

    worse_arm = np.maximum(
        single["normal_max_log"].to_numpy(dtype=np.float32)[ia],
        single["normal_max_log"].to_numpy(dtype=np.float32)[ib],
    )
    feats["and_gate_gain"] = worse_arm - feats["normal_and_max_log"].to_numpy()

    # ---------- tumour co-expression (cancer-type averages) ----------
    if verbose:
        print("  tumour co-expression (cancer-type means)")
    tum_log = tumour_profiles(universe)
    TL = tum_log.to_numpy(dtype=np.float32)
    TRAW = np.expm1(TL)
    feats["tumour_and_max_log"] = np.log1p(P.max_of_min(TRAW, ia, ib))
    feats["tumour_double_pos_n"] = P.both_above(TRAW, ia, ib, C.EXPR_EXPRESSED)
    feats["tumour_profile_corr"] = P.row_corr(TL, ia, ib)

    # ---------- tumour co-expression (individual cell lines) ----------
    # Stricter than the above: asks whether individual tumour samples carry
    # both antigens, not whether the cancer-type averages do.
    if verbose:
        print("  tumour co-expression (1,206 individual cell lines)")
    lines = D.load_tumour_lines().reindex(genes).fillna(0.0)
    L = np.ascontiguousarray(lines.to_numpy(dtype=np.float32))
    n_lines = max(L.shape[1], 1)
    feats["line_double_pos_frac"] = P.both_above(L, ia, ib, C.EXPR_EXPRESSED) / n_lines
    feats["line_either_pos_frac"] = P.either_above(L, ia, ib, C.EXPR_EXPRESSED) / n_lines
    feats["line_and_max_log"] = np.log1p(P.max_of_min(L, ia, ib))
    feats["line_profile_corr"] = P.row_corr(np.log1p(L), ia, ib)
    del L, lines

    # ---------- pair-level therapeutic window ----------
    feats["pair_window_log"] = feats["tumour_and_max_log"] - feats["normal_and_max_log"]
    feats["pair_window_vital_log"] = (
        feats["tumour_and_max_log"] - feats.get("normal_and_vital_max_log", 0.0)
    )

    # ---------- functional overlap ----------
    if verbose:
        print("  functional term overlap")
    annot = D.load_annotations().set_index("gene").reindex(genes)
    for src, dst in [("Biological process", "bp_jaccard"), ("Molecular function", "mf_jaccard")]:
        sets = annot[src].map(D.split_terms).to_list()
        vals = np.full(len(ia), np.nan, dtype=np.float32)
        for k, (i, j) in enumerate(zip(ia, ib)):
            a, b = sets[i], sets[j]
            if a and b:
                union = len(a | b)
                vals[k] = len(a & b) / union if union else np.nan
        feats[dst] = vals

    # ---------- tumour-arm-resolved features ----------
    role = universe.set_index("gene")["role"].reindex(genes).to_numpy()
    if verbose:
        print("  resolving features to the tumour arm")
    arm = _taa_arm_features(single, ia, ib, role)
    for c in arm.columns:
        feats[c] = arm[c].to_numpy()

    # ---------- role combination ----------
    combo = np.array(["|".join(sorted((role[i], role[j]))) for i, j in zip(ia, ib)])
    feats["role_combo"] = combo
    for lab in ["TAA|effector", "TAA|checkpoint", "TAA|TAA",
                "checkpoint|checkpoint", "checkpoint|effector"]:
        feats[f"role__{lab.replace('|', '_')}"] = (combo == lab).astype(np.float32)

    # ---------- single-cell co-detection ----------
    # Defragment: the frame was built by ~200 successive column assignments.
    feats = feats.copy().set_index("pair")
    if sc_features is not None and len(sc_features):
        if verbose:
            print("  joining single-cell co-detection features")
        feats = feats.join(sc_features, how="left")
    return feats.reset_index()


def label_pairs(pairs: pd.DataFrame) -> pd.Series:
    """1 for documented bispecific pairs, 0 for everything else.

    Note these zeros are *unlabelled*, not verified negatives: the great
    majority of pairs have simply never been tried. Downstream evaluation
    treats this as positive-unlabelled learning.
    """
    from .knowledge import known_pairs_df

    known = set(known_pairs_df().pair)
    return pairs.pair.isin(known).astype(int)
