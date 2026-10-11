"""End-to-end pipeline stages, in dependency order.

Each stage caches its output under data/processed/ so a rerun resumes rather
than recomputing. Run via `python run.py` at the project root.
"""

from __future__ import annotations

import time

import pandas as pd

from . import config as C
from . import features_pair, features_sc, features_single, score, targets
from .knowledge import known_pairs_df


def _stamp(msg: str, t0: float | None = None) -> float:
    now = time.time()
    if t0 is None:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}")
    else:
        print(f"[{time.strftime('%H:%M:%S')}] {msg} ({now - t0:.1f}s)")
    return now


def stage_universe(force: bool = False) -> pd.DataFrame:
    path = C.DATA_PROC / "universe.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    t0 = _stamp("building target universe")
    uni = targets.build_universe()
    uni.to_parquet(path)
    uni.to_csv(C.DATA_PROC / "universe.csv", index=False)
    targets.benchmark_priors(uni).to_csv(C.RESULTS / "prior_benchmark.csv", index=False)
    audit = targets.audit_surface_filter(verbose=False)
    pd.DataFrame([{k: (", ".join(v) if isinstance(v, list) else v)
                   for k, v in audit.items()}]).to_csv(
        C.RESULTS / "surface_filter_audit.csv", index=False)
    print(f"  surface filter on {audit['known_genes']} validated targets: "
          f"{audit['high']} high, {audit['medium']} medium, "
          f"{audit['rejected']} rejected ({', '.join(audit['rejected_genes'])})")
    _stamp("universe done", t0)
    return uni


def stage_pairs(universe: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    path = C.DATA_PROC / "pairs.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    t0 = _stamp("building pair space")
    pairs = targets.build_pair_space(universe)
    pairs.to_parquet(path)
    _stamp("pair space done", t0)
    return pairs


def stage_single(universe: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    path = C.DATA_PROC / "single_features.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    t0 = _stamp("building single-arm features")
    single = features_single.build(universe)
    single.to_parquet(path)
    _stamp(f"single-arm features done: {single.shape}", t0)
    return single


def stage_sc(universe: pd.DataFrame, pairs: pd.DataFrame, force: bool = False) -> pd.DataFrame:
    path = C.DATA_PROC / "sc_features.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    t0 = _stamp("single-cell co-detection across tumour datasets")
    stats = features_sc.collect(list(universe.gene))
    if not stats:
        print("  no usable single-cell datasets found")
        empty = pd.DataFrame(index=pd.Index(pairs.pair, name="pair"))
        empty.to_parquet(path)
        return empty
    sc = features_sc.pair_sc_features(stats, pairs)
    summary = pd.DataFrame([
        {"dataset": s["name"], "disease": s["disease"],
         "malignant_cells": s["malignant"]["n_cells"],
         "other_cells": s.get("stromal", {}).get("n_cells", 0),
         "universe_genes_present": len(s["genes"])}
        for s in stats
    ])
    summary.to_csv(C.RESULTS / "scrna_datasets.csv", index=False)
    print(summary.to_string(index=False))
    sc.to_parquet(path)
    _stamp(f"single-cell features done: {sc.shape}", t0)
    return sc


def stage_features(universe, pairs, single, sc, force: bool = False) -> pd.DataFrame:
    path = C.DATA_PROC / "pair_features.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    t0 = _stamp("assembling pair features")
    feats = features_pair.build(universe, pairs, single, sc)
    feats.to_parquet(path)
    _stamp(f"pair features done: {feats.shape}", t0)
    return feats


def stage_score(features: pd.DataFrame, pairs: pd.DataFrame, force: bool = False):
    path = C.DATA_PROC / "scored.parquet"
    t0 = _stamp("scoring pairs with the role-aware rubric")
    scored = score.score_pairs(features)
    labels = features_pair.label_pairs(features[["pair"]])
    validation = score.validate_against_known(scored, labels)
    scored.to_parquet(path)
    validation.to_csv(C.RESULTS / "rubric_validation.csv", index=False)
    score.component_table().to_csv(C.RESULTS / "rubric_components.csv", index=False)
    print(validation.to_string(index=False))
    _stamp("scoring done", t0)
    return scored, validation


def run_all(force: bool = False):
    universe = stage_universe(force)
    pairs = stage_pairs(universe, force)
    single = stage_single(universe, force)
    sc = stage_sc(universe, pairs, force)
    features = stage_features(universe, pairs, single, sc, force)
    scored, validation = stage_score(features, pairs, force)
    return {
        "universe": universe,
        "pairs": pairs,
        "single": single,
        "sc": sc,
        "features": features,
        "scored": scored,
        "validation": validation,
        "known": known_pairs_df(),
    }
