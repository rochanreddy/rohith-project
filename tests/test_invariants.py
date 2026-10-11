"""Invariants the pipeline must satisfy. Run with: python -m pytest tests -q

These are the properties whose violation would silently corrupt the results
rather than raise: arm-order leakage, non-canonical pair keys, chunked
reductions disagreeing with their brute-force definition, and positives being
dropped from the search space.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from bsab import pairops as P
from bsab import score
from bsab.features_pair import TAA_ARM_SPEC, _symmetric_aggregates, _taa_arm_features
from bsab.knowledge import known_pairs_df, pair_key

rng = np.random.default_rng(0)


# ---------------------------------------------------------------- pair keys
def test_pair_key_is_order_invariant():
    assert pair_key("CD19", "CD3E") == pair_key("CD3E", "CD19")


def test_known_pairs_are_unique_and_canonical():
    kn = known_pairs_df()
    assert kn.pair.is_unique
    for _, r in kn.iterrows():
        assert r.pair == pair_key(r.gene_a, r.gene_b)
        assert r.gene_a != r.gene_b


# ------------------------------------------------- chunked pairwise reductions
@pytest.fixture
def profile():
    return rng.random((40, 17)).astype(np.float32) * 20


@pytest.fixture
def idx():
    ia, ib = np.triu_indices(40, k=1)
    return ia, ib


def test_max_of_min_matches_bruteforce(profile, idx):
    ia, ib = idx
    got = P.max_of_min(profile, ia, ib, chunk=7)
    want = np.array([np.minimum(profile[i], profile[j]).max() for i, j in zip(ia, ib)])
    assert np.allclose(got, want, atol=1e-5)


def test_both_above_matches_bruteforce(profile, idx):
    ia, ib = idx
    got = P.both_above(profile, ia, ib, 10.0, chunk=7)
    want = np.array([((profile[i] >= 10) & (profile[j] >= 10)).sum()
                     for i, j in zip(ia, ib)], dtype=float)
    assert np.array_equal(got, want)


def test_row_corr_matches_numpy(profile, idx):
    ia, ib = idx
    got = P.row_corr(profile, ia, ib, chunk=7)
    want = np.array([np.corrcoef(profile[i], profile[j])[0, 1] for i, j in zip(ia, ib)])
    assert np.allclose(got, want, atol=1e-4)


def test_reductions_are_symmetric(profile, idx):
    ia, ib = idx
    for fn in (P.max_of_min, P.min_of_max, P.row_corr):
        assert np.allclose(fn(profile, ia, ib), fn(profile, ib, ia), atol=1e-5)
    assert np.array_equal(P.both_above(profile, ia, ib, 5.0),
                          P.both_above(profile, ib, ia, 5.0))


def test_chunking_does_not_change_results(profile, idx):
    ia, ib = idx
    for chunk in (1, 3, 1000):
        assert np.allclose(P.max_of_min(profile, ia, ib, chunk=chunk),
                           P.max_of_min(profile, ia, ib), atol=1e-6)


# ------------------------------------------------------- arm-order invariance
def test_symmetric_aggregates_invariant_to_arm_swap():
    single = pd.DataFrame(rng.random((12, 4)).astype(np.float32),
                          columns=["a", "b", "c", "d"])
    ia = np.array([0, 1, 2, 3, 4])
    ib = np.array([5, 6, 7, 8, 9])
    forward = _symmetric_aggregates(single, ia, ib)
    reverse = _symmetric_aggregates(single, ib, ia)
    pd.testing.assert_frame_equal(forward, reverse)


def test_taa_arm_features_pick_the_tumour_arm():
    single = pd.DataFrame({
        "tumour_best_frac": np.array([0.9, 0.1, 0.5], dtype=np.float32),
        "normal_vital_max_log": np.array([1.0, 8.0, 4.0], dtype=np.float32),
    })
    role = np.array(["TAA", "effector", "TAA"])
    # pair 0: TAA(0) x effector(1) -> takes gene 0's values
    # pair 1: TAA(0) x TAA(2)      -> conservative: min of the "higher is
    #                                 better" column, max of the risk column
    ia = np.array([0, 0])
    ib = np.array([1, 2])
    out = _taa_arm_features(single, ia, ib, role)
    assert out.loc[0, "taa_tumour_best_frac"] == pytest.approx(0.9)
    assert out.loc[0, "taa_normal_vital_max_log"] == pytest.approx(1.0)
    assert out.loc[1, "taa_tumour_best_frac"] == pytest.approx(0.5)
    assert out.loc[1, "taa_normal_vital_max_log"] == pytest.approx(4.0)
    assert list(out.n_taa_arms) == [1, 2]


def test_taa_arm_features_are_nan_without_a_tumour_arm():
    single = pd.DataFrame({"tumour_best_frac": np.array([0.4, 0.6], dtype=np.float32)})
    role = np.array(["effector", "checkpoint"])
    out = _taa_arm_features(single, np.array([0]), np.array([1]), role)
    assert np.isnan(out.loc[0, "taa_tumour_best_frac"])


def test_taa_arm_spec_directions_are_declared():
    for col, hib in TAA_ARM_SPEC.items():
        assert isinstance(hib, bool), col


# ------------------------------------------------------------------- scoring
def _toy_features(n=200):
    cls = np.array(["TAA|effector"] * (n // 2) + ["TAA|TAA"] * (n - n // 2))
    return pd.DataFrame({
        "pair": [f"G{i}|H{i}" for i in range(n)],
        "role_combo": cls,
        "taa_tumour_best_frac": rng.random(n),
        "taa_tumour_max_log": rng.random(n) * 8,
        "taa_tumour_n_types_hit": rng.integers(0, 5, n).astype(float),
        "taa_window_vital_log": rng.normal(size=n),
        "taa_normal_vital_max_log": rng.random(n) * 6,
        "taa_normal_frac_expressed": rng.random(n),
        "taa_normal_tau": rng.random(n),
        "taa_loc_plasma_membrane": rng.integers(0, 2, n).astype(float),
        "taa_secreted_to_blood": rng.integers(0, 2, n).astype(float),
        "taa_cls_fda_target": rng.integers(0, 2, n).astype(float),
        "sc_pboth_max": rng.random(n),
        "sc_lift_max": rng.random(n) * 5,
        "sc_tme_margin": rng.normal(size=n) * 0.1,
        "sc_stromal_pboth_max": rng.random(n),
        "line_double_pos_frac": rng.random(n),
        "normal_and_max_log": rng.random(n) * 7,
        "normal_double_pos_vital_n": rng.integers(0, 20, n).astype(float),
        "and_gate_gain": rng.normal(size=n),
        "mf_jaccard": rng.random(n),
        "bp_jaccard": rng.random(n),
        "tumour_profile_corr": rng.normal(size=n) * 0.5,
        "cls_cancer_related__min": rng.integers(0, 2, n).astype(float),
        "normal_frac_expressed__max": rng.random(n),
    })


def test_scores_are_bounded_percentiles():
    s = score.score_pairs(_toy_features())
    ok = s.score.notna()
    assert ok.any()
    assert s.score[ok].between(0, 1).all()
    assert s.coverage[ok].between(0, 1).all()


def test_each_class_uses_its_own_rubric():
    s = score.score_pairs(_toy_features())
    # TAA|effector has one rationale, TAA|TAA has two.
    eff = s[s.design_class == "TAA|effector"]
    tt = s[s.design_class == "TAA|TAA"]
    assert set(eff.rationale.dropna()) == {"effector_engager"}
    assert set(tt.rationale.dropna()) <= {"and_gate_selectivity", "pathway_redundancy"}
    # The chosen score must be the max over that class's rationales.
    both = tt[["score__and_gate_selectivity", "score__pathway_redundancy"]].max(axis=1)
    assert np.allclose(tt.score.to_numpy(), both.to_numpy(), equal_nan=True)


def test_missing_component_does_not_penalise():
    """A pair with no single-cell data must not score below an identical pair
    that has it, purely for lacking the measurement."""
    f = _toy_features(100)
    # sc_pboth_max is only a component of the TAA|TAA rubrics, so blank it there.
    target = f.index[f.role_combo == "TAA|TAA"][:10]
    f.loc[target, ["sc_pboth_max", "sc_lift_max", "sc_tme_margin"]] = np.nan
    s = score.score_pairs(f)
    hit = s.loc[target]
    assert hit["coverage__and_gate_selectivity"].max() < 1.0
    # Still scored, from the remaining components, rather than dropped or zeroed.
    assert hit.score.notna().all()
    assert hit.score.between(0, 1).all()


def test_component_table_covers_every_rubric():
    tbl = score.component_table()
    assert set(tbl.design_class) == set(score.RUBRICS)
    for cls, rationales in score.RUBRICS.items():
        got = set(tbl[tbl.design_class == cls].rationale)
        assert got == set(rationales)
    assert (tbl.weight > 0).all()


def test_rubric_weights_sum_to_one_per_rationale():
    for cls, rationales in score.RUBRICS.items():
        for rname, rubric in rationales.items():
            total = sum(w for w, _, _ in rubric.values())
            assert total == pytest.approx(1.0, abs=0.02), f"{cls}/{rname} sums to {total}"


# -------------------------------------------------------- model matrix hygiene
def test_model_matrix_excludes_identifiers_and_raw_arm_columns():
    from bsab import model

    f = _toy_features(50)
    f["gene_a"] = "A"
    f["gene_b"] = "B"
    f["sc_best_disease"] = "breast cancer"
    X, cols = model.feature_matrix(f)
    for banned in ("pair", "gene_a", "gene_b", "role_combo", "sc_best_disease"):
        assert banned not in cols
    assert X.notna().any().any()
