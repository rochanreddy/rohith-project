"""Supervised pairwise model, with evaluation designed to expose its limits.

The honest problem statement: 64 documented bispecific pairs, 48 of them
sharing CD3E as one arm, against 232,450 unlabelled pairs that are *not*
verified negatives. Three things follow, and each is handled explicitly:

1. **Arm leakage.** Random k-fold would put CD19|CD3E in train and CD22|CD3E in
   test, and a model scoring "contains CD3E" would look excellent. All
   evaluation is therefore grouped by the shared arm: every pair involving a
   given gene sits in exactly one fold (`GroupKFold` on the rarer arm).

2. **CD3E dominance.** Because T-cell engagers dominate the literature, we also
   evaluate on the 16 non-CD3E pairs only, after training with every CD3E pair
   removed. That is the question a discovery team actually cares about: can the
   model find a good pair it was not shown an example of?

3. **Unlabelled negatives.** Reported as enrichment over chance and recall@k,
   not as accuracy. ROC-AUC against unlabelled negatives is reported but should
   be read as an upper bound: some "negatives" are simply untried pairs.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C

DROP_COLS = {"pair", "gene_a", "gene_b", "block", "role_combo", "sc_best_disease"}


def feature_matrix(features: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    cols = [c for c in features.columns
            if c not in DROP_COLS and pd.api.types.is_numeric_dtype(features[c])]
    X = features[cols].astype("float32")
    return X, cols


def _group_labels(pairs: pd.DataFrame, features_index) -> np.ndarray:
    """Group each pair by its rarer arm, so no gene spans two folds.

    Using the rarer arm keeps the hub genes (CD3E and friends) from collapsing
    everything into one group.
    """
    counts = pd.concat([pairs.gene_a, pairs.gene_b]).value_counts()
    ga = pairs.gene_a.to_numpy()
    gb = pairs.gene_b.to_numpy()
    ca = counts.reindex(ga).to_numpy()
    cb = counts.reindex(gb).to_numpy()
    return np.where(ca <= cb, ga, gb)


def _lgbm(seed: int = C.RANDOM_SEED):
    import lightgbm as lgb

    return lgb.LGBMClassifier(
        n_estimators=400,
        learning_rate=0.05,
        num_leaves=31,
        min_child_samples=30,
        subsample=0.8,
        subsample_freq=1,
        colsample_bytree=0.6,
        reg_lambda=5.0,
        is_unbalance=True,
        random_state=seed,
        n_jobs=-1,
        verbose=-1,
    )


def _metrics(y: np.ndarray, s: np.ndarray) -> dict:
    from sklearn.metrics import average_precision_score, roc_auc_score

    npos, n = int(y.sum()), len(y)
    res: dict[str, object] = {"n": n, "n_pos": npos}
    if npos == 0 or npos == n:
        return res
    res["roc_auc"] = round(float(roc_auc_score(y, s)), 4)
    res["average_precision"] = round(float(average_precision_score(y, s)), 4)
    order = np.argsort(-s)
    ys = y[order]
    for k in (100, 500, 1000, 5000):
        if k > n:
            continue
        hit = int(ys[:k].sum())
        res[f"recall@{k}"] = round(hit / npos, 3)
        res[f"enrichment@{k}"] = round((hit / k) / (npos / n), 2)
    return res


def cross_validate(features: pd.DataFrame, pairs: pd.DataFrame, labels: pd.Series,
                   n_splits: int = 5, exclude_cd3e: bool = False,
                   verbose: bool = True) -> dict:
    """Arm-grouped cross-validation over the full pair space.

    If `exclude_cd3e`, every pair involving CD3E is dropped from *training*
    folds, and evaluation is restricted to non-CD3E pairs - the generalisation
    test that matters.
    """
    from sklearn.model_selection import GroupKFold

    pairs = pairs.set_index("pair").loc[features["pair"]].reset_index()
    X, cols = feature_matrix(features)
    y = labels.to_numpy().astype(int)
    groups = _group_labels(pairs, features.index)

    has_cd3e = ((pairs.gene_a == "CD3E") | (pairs.gene_b == "CD3E")).to_numpy()

    oof = np.full(len(y), np.nan)
    gkf = GroupKFold(n_splits=n_splits)
    for fold, (tr, te) in enumerate(gkf.split(X, y, groups)):
        if exclude_cd3e:
            tr = tr[~has_cd3e[tr]]
        if y[tr].sum() == 0:
            continue
        clf = _lgbm()
        clf.fit(X.iloc[tr], y[tr])
        oof[te] = clf.predict_proba(X.iloc[te])[:, 1]
        if verbose:
            print(f"    fold {fold}: train {len(tr):,} ({int(y[tr].sum())} pos) "
                  f"-> test {len(te):,} ({int(y[te].sum())} pos)")

    mask = ~np.isnan(oof)
    if exclude_cd3e:
        mask &= ~has_cd3e
    out = _metrics(y[mask], oof[mask])
    out["setting"] = "arm-grouped CV, CD3E removed from training" if exclude_cd3e else "arm-grouped CV"
    return {"metrics": out, "oof": oof, "mask": mask, "columns": cols}


def random_cv_for_contrast(features: pd.DataFrame, labels: pd.Series,
                           n_splits: int = 5) -> dict:
    """Deliberately naive random k-fold, to quantify how much arm leakage inflates scores.

    Included because this is the number a careless version of this project
    would report.
    """
    from sklearn.model_selection import StratifiedKFold

    X, _ = feature_matrix(features)
    y = labels.to_numpy().astype(int)
    oof = np.full(len(y), np.nan)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=C.RANDOM_SEED)
    for tr, te in skf.split(X, y):
        clf = _lgbm()
        clf.fit(X.iloc[tr], y[tr])
        oof[te] = clf.predict_proba(X.iloc[te])[:, 1]
    res = _metrics(y, oof)
    res["setting"] = "random k-fold (leaky - for contrast only)"
    return {"metrics": res, "oof": oof}


def fit_full(features: pd.DataFrame, labels: pd.Series):
    X, cols = feature_matrix(features)
    clf = _lgbm()
    clf.fit(X, labels.to_numpy().astype(int))
    imp = (pd.DataFrame({"feature": cols, "gain": clf.booster_.feature_importance("gain")})
           .sort_values("gain", ascending=False).reset_index(drop=True))
    imp["gain_frac"] = imp.gain / imp.gain.sum()
    return clf, imp, cols


def contributions(clf, features: pd.DataFrame, cols: list[str],
                  rows: np.ndarray | None = None) -> pd.DataFrame:
    """Per-pair additive feature contributions (LightGBM SHAP) for explanations."""
    X = features[cols].astype("float32")
    if rows is not None:
        X = X.iloc[rows]
    raw = clf.booster_.predict(X, pred_contrib=True)
    contrib = pd.DataFrame(raw[:, :-1], columns=cols, index=X.index)
    contrib["_bias"] = raw[:, -1]
    return contrib
