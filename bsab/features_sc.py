"""Single-cell co-expression features from real tumour scRNA-seq datasets.

Why this matters for bispecifics: a pair can look good on bulk or averaged data
while the two antigens actually sit on *different* cells within the tumour. An
avidity-driven or AND-gated bispecific only works when both antigens are on the
*same* cell, so co-detection has to be measured cell by cell. This is the one
part of the pipeline that needs no labels at all.

Detection is defined as value > 0, which is invariant to per-cell library-size
normalisation and to log1p, so it is comparable across datasets processed
differently. Datasets whose matrix has been z-scored (negative values) fall
back to the stored raw counts.

Caveat carried into the report: scRNA-seq dropout means a zero is not proof of
absence, so co-detection rates are systematically *under*-estimates. They are
useful for ranking pairs against each other, not as absolute surface-protein
co-expression, and mRNA is a proxy for protein in any case.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import config as C

MALIGNANT_LABELS = {"malignant cell", "neoplastic cell"}


def _detection_matrix(adata, genes: list[str]):
    """Binary cells x genes detection matrix restricted to `genes`."""
    source, X = adata, adata.X
    probe = X[:50]
    mn = probe.min() if sp.issparse(probe) else np.nanmin(probe)
    if mn < 0 and adata.raw is not None:
        source, X = adata.raw, adata.raw.X

    if "feature_name" in source.var:
        names = source.var["feature_name"].astype(str).to_numpy()
    else:
        names = source.var.index.astype(str).to_numpy()

    pos: dict[str, int] = {}
    for i, g in enumerate(names):
        pos.setdefault(g, i)

    present = [g for g in genes if g in pos]
    idx = [pos[g] for g in present]
    sub = X[:, idx]
    if not sp.issparse(sub):
        sub = sp.csr_matrix(sub)
    return sp.csr_matrix((sub > 0).astype(np.float32)), present


def _codetect(B: sp.csr_matrix) -> dict:
    """Pairwise co-detection counts for a binary cells x genes matrix."""
    both = (B.T @ B)
    both = np.asarray(both.todense(), dtype=np.float32) if sp.issparse(both) else np.asarray(both, dtype=np.float32)
    return {"n_cells": int(B.shape[0]), "both": both, "single": np.diag(both).copy()}


def dataset_stats(path, genes: list[str], verbose: bool = True) -> dict | None:
    """Co-detection stats for malignant and non-malignant cells of one dataset."""
    import anndata as ad

    adata = ad.read_h5ad(path)
    labels = adata.obs["cell_type"].astype(str)
    is_mal = labels.isin(MALIGNANT_LABELS).to_numpy()
    if is_mal.sum() < C.SC_MIN_MALIGNANT_CELLS:
        if verbose:
            print(f"  {path.stem}: only {int(is_mal.sum())} malignant cells, skipped")
        return None

    B, present = _detection_matrix(adata, genes)
    disease = str(labels.index.to_series().pipe(lambda _: adata.obs["disease"].astype(str).mode().iloc[0]))

    out = {
        "name": path.stem,
        "disease": disease,
        "genes": present,
        "gidx": {g: i for i, g in enumerate(present)},
        "malignant": _codetect(sp.csr_matrix(B[is_mal])),
    }
    if (~is_mal).sum() >= 100:
        out["stromal"] = _codetect(sp.csr_matrix(B[~is_mal]))
        out["stromal_composition"] = labels[~is_mal].value_counts().head(8).to_dict()
    if verbose:
        print(f"  {path.stem}: {int(is_mal.sum())} malignant / {int((~is_mal).sum())} other cells, "
              f"{len(present)} of {len(genes)} universe genes present ({disease})")
    del adata
    return out


def collect(genes: list[str], paths=None, verbose: bool = True) -> list[dict]:
    paths = list(paths) if paths else sorted(C.SCRNA_DIR.glob("*.h5ad"))
    stats = []
    for p in paths:
        s = dataset_stats(p, genes, verbose=verbose)
        if s is not None:
            stats.append(s)
    return stats


def pair_sc_features(stats: list[dict], pairs: pd.DataFrame) -> pd.DataFrame:
    """Per-pair single-cell co-expression features, aggregated over datasets.

    Per dataset, over malignant cells only:
      p_a, p_b   detection rate of each arm
      p_both     fraction of malignant cells detecting both
      jaccard    p_both / (p_a + p_b - p_both)
      lift       p_both / (p_a * p_b) - above 1 means the two antigens coincide
                 on the same cells more than independent expression would give
    and over the non-malignant cells of the same tumour:
      stromal_p_both  double-positive rate in the microenvironment, a direct
                      in-situ read-out of on-target/off-tumour risk
    """
    n = len(pairs)
    ga = pairs.gene_a.to_numpy()
    gb = pairs.gene_b.to_numpy()

    acc = {k: np.full(n, np.nan, dtype=np.float32) for k in
           ["pboth_max", "pboth_mean", "jaccard_max", "lift_max", "lift_mean",
            "stromal_pboth_max", "tme_margin", "n_datasets"]}
    best_disease = np.array([None] * n, dtype=object)
    per_dataset: dict[str, np.ndarray] = {}

    pboth_stack, jac_stack, lift_stack, strom_stack = [], [], [], []

    for s in stats:
        gidx = s["gidx"]
        ja = np.array([gidx.get(g, -1) for g in ga])
        jb = np.array([gidx.get(g, -1) for g in gb])
        ok = (ja >= 0) & (jb >= 0)

        mal = s["malignant"]
        nc = mal["n_cells"]
        pa = np.full(n, np.nan, dtype=np.float32)
        pb = np.full(n, np.nan, dtype=np.float32)
        pboth = np.full(n, np.nan, dtype=np.float32)
        pa[ok] = mal["single"][ja[ok]] / nc
        pb[ok] = mal["single"][jb[ok]] / nc
        pboth[ok] = mal["both"][ja[ok], jb[ok]] / nc

        union = pa + pb - pboth
        with np.errstate(invalid="ignore", divide="ignore"):
            jac = np.where(union > 0, pboth / union, 0.0)
            lift = np.where((pa * pb) > 0, pboth / (pa * pb), 0.0)

        per_dataset[f"sc_pboth__{s['name']}"] = pboth
        pboth_stack.append(pboth)
        jac_stack.append(jac.astype(np.float32))
        lift_stack.append(lift.astype(np.float32))

        if "stromal" in s:
            st = s["stromal"]
            sp_both = np.full(n, np.nan, dtype=np.float32)
            sp_both[ok] = st["both"][ja[ok], jb[ok]] / st["n_cells"]
            strom_stack.append(sp_both)

    if not pboth_stack:
        return pd.DataFrame(index=pd.Index(pairs.pair, name="pair"))

    PB = np.vstack(pboth_stack)
    JA = np.vstack(jac_stack)
    LI = np.vstack(lift_stack)

    with np.errstate(invalid="ignore"):
        acc["n_datasets"] = np.sum(~np.isnan(PB), axis=0).astype(np.float32)
        acc["pboth_max"] = np.nanmax(PB, axis=0)
        acc["pboth_mean"] = np.nanmean(PB, axis=0)
        acc["jaccard_max"] = np.nanmax(JA, axis=0)
        acc["lift_max"] = np.nanmax(LI, axis=0)
        acc["lift_mean"] = np.nanmean(LI, axis=0)
        if strom_stack:
            SB = np.vstack(strom_stack)
            acc["stromal_pboth_max"] = np.nanmax(SB, axis=0)
            acc["tme_margin"] = acc["pboth_max"] - acc["stromal_pboth_max"]

    # Which tumour type gives the strongest same-cell co-detection.
    diseases = np.array([s["disease"] for s in stats], dtype=object)
    with np.errstate(invalid="ignore"):
        allnan = np.all(np.isnan(PB), axis=0)
        arg = np.zeros(PB.shape[1], dtype=int)
        arg[~allnan] = np.nanargmax(PB[:, ~allnan], axis=0)
    best_disease = np.where(allnan, None, diseases[arg])

    out = pd.DataFrame({f"sc_{k}": v for k, v in acc.items()})
    for k, v in per_dataset.items():
        out[k] = v
    out["sc_best_disease"] = best_disease
    out.index = pd.Index(pairs.pair, name="pair")
    return out
