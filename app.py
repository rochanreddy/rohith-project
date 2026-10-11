#!/usr/bin/env python
"""Browsable UI over the bispecific target-pair results.

    .venv/Scripts/streamlit run app.py        # -> http://localhost:8501

Reads whatever `run.py` has produced in results/. Nothing is computed here
except the per-pair explanation text, which is generated on demand from the
cached feature table so it always matches the numbers on screen.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"
PROC = ROOT / "data" / "processed"

st.set_page_config(page_title="Bispecific target pairs", page_icon="🧬", layout="wide")

CAVEAT = (
    "Target triage to prioritise experiments — **not** evidence of efficacy or safety. "
    "Expression is mRNA, not surface protein; scRNA-seq dropout makes co-detection a "
    "systematic under-estimate; cell lines lack microenvironment. Nothing here should "
    "inform a clinical decision."
)


# ----------------------------------------------------------------- loading
@st.cache_data(show_spinner=False)
def load_csv(name: str) -> pd.DataFrame:
    p = RESULTS / name
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


@st.cache_data(show_spinner=False)
def load_ranked() -> pd.DataFrame:
    return load_csv("ranked_pairs.csv")


@st.cache_data(show_spinner=False)
def load_features() -> pd.DataFrame:
    p = PROC / "pair_features.parquet"
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()


@st.cache_data(show_spinner=False)
def load_universe() -> pd.DataFrame:
    p = PROC / "universe.parquet"
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()


@st.cache_data(show_spinner=False)
def load_normal() -> pd.DataFrame:
    p = PROC / "normal_sc_wide.parquet"
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()


@st.cache_data(show_spinner=False)
def load_scored() -> pd.DataFrame:
    """Per-component percentiles. These live only in the parquet: adding ~45
    float columns to ranked_pairs.csv would push it past 80 MB."""
    p = PROC / "scored.parquet"
    return pd.read_parquet(p) if p.exists() else pd.DataFrame()


ranked = load_ranked()
if ranked.empty:
    st.title("Bispecific antibody target pairs")
    st.error("No results found. Run the pipeline first:  `python run.py`")
    st.stop()

# ----------------------------------------------------------------- sidebar
st.sidebar.title("Filters")
classes = sorted(ranked.design_class.dropna().unique())
pick_class = st.sidebar.multiselect("Design class", classes, default=classes)

rationales = sorted(ranked.rationale.dropna().unique()) if "rationale" in ranked else []
pick_rationale = st.sidebar.multiselect("Rationale", rationales, default=rationales)

show_what = st.sidebar.radio("Show", ["Candidates only", "Documented only", "Everything"])
min_score = st.sidebar.slider("Minimum score", 0.0, 1.0, 0.0, 0.01)
min_cov = st.sidebar.slider("Minimum evidence coverage", 0.0, 1.0, 0.0, 0.05,
                            help="Fraction of the rubric's weight that had data. "
                                 "Low coverage means the score rests on few components.")
gene_query = st.sidebar.text_input("Gene contains", "", placeholder="e.g. CD3E, ERBB2")
hide_low_conf = st.sidebar.checkbox("Hide low-confidence classes", value=False,
                                    help="checkpoint x checkpoint and checkpoint x effector: "
                                         "expression data says very little about these.")
surface_only = st.sidebar.checkbox(
    "Well-evidenced surface targets only", value=True,
    help="HPA localisation cannot tell which FACE of the membrane a protein sits on, so "
         "cytoplasmic adaptors (PPP1R18, PLEKHO2) are admitted as plausible targets and can "
         "reach rank 1. Leaving this on keeps the tumour arm to proteins with positive "
         "surface/extracellular evidence. Turn it off to see everything, including artifacts.")

view = ranked[ranked.design_class.isin(pick_class)].copy()
if pick_rationale and "rationale" in view:
    view = view[view.rationale.isin(pick_rationale)]
if show_what == "Candidates only":
    view = view[~view.is_documented.astype(bool)]
elif show_what == "Documented only":
    view = view[view.is_documented.astype(bool)]
view = view[view.score.fillna(0) >= min_score]
if "coverage" in view:
    view = view[view.coverage.fillna(0) >= min_cov]
if gene_query.strip():
    view = view[view.pair.str.contains(gene_query.strip().upper(), na=False)]
if hide_low_conf and "low_confidence" in view:
    view = view[~view.low_confidence.astype(bool)]
if surface_only and "taa_surface_confidence" in view:
    view = view[view.taa_surface_confidence.fillna(0) >= 1.0]

# ----------------------------------------------------------------- header
st.title("Bispecific antibody target pairs")
st.caption(CAVEAT)

c1, c2, c3, c4 = st.columns(4)
c1.metric("Pairs scored", f"{len(ranked):,}")
c2.metric("Matching filters", f"{len(view):,}")
c3.metric("Documented pairs", int(ranked.is_documented.astype(bool).sum()))
c4.metric("Design classes", ranked.design_class.nunique())

tabs = st.tabs(["⭐ Start here", "Ranked pairs", "Pair detail", "By indication",
                "Validation", "Figures", "Report"])

# ----------------------------------------------------------------- story
with tabs[0]:
    st.markdown("### What problem is this solving?")
    st.markdown(
        "A **bispecific antibody** is a drug with two arms. One arm grabs a **cancer cell**, "
        "the other grabs a **T-cell** — physically dragging the immune system onto the tumour. "
        "Several are approved drugs today.\n\n"
        "The catch: about **3,000 human proteins** sit on the cell surface where an antibody "
        "can reach them, which makes **millions of possible two-target combinations**. "
        "Which pair deserves years of laboratory work?\n\n"
        "Two requirements decide it, and they pull against each other:"
    )
    a, b = st.columns(2)
    a.success("**It works** if both targets are on the tumour — ideally on the *same* tumour cell.")
    b.error("**It is unsafe** if healthy tissue also carries both, because the drug attacks that too.")

    for img, cap in [
        ("story_pipeline.png", None),
        ("story_same_cell.png",
         "Averaged expression data cannot tell these two tumours apart — only "
         "cell-by-cell measurement can, and it decides whether the antibody can physically work."),
    ]:
        p = FIGURES / img
        if p.exists():
            st.image(str(p), width='stretch')
            if cap:
                st.caption(cap)

    st.markdown("---")
    st.markdown("### Does it actually find real drugs?")
    st.markdown(
        "The honest test: hide nothing, and check whether the method rediscovers bispecifics "
        "that already exist. **BCMA** is the clearest case — the target of *teclistamab*, an "
        "approved myeloma drug."
    )
    m1, m2, m3 = st.columns(3)
    m1.metric("BCMA in myeloma", "rank 1", "found without being told")
    m2.metric("Present in myeloma samples", "100%", "vs ~0% in other cancers")
    m3.metric("Healthy tissue carrying it", "plasma cells", "a lineage patients survive losing")
    for img in ("story_example_tumour.png", "story_example_normal.png"):
        p = FIGURES / img
        if p.exists():
            st.image(str(p), width='stretch')
    st.markdown(
        "That second chart is the whole argument. BCMA *does* have high expression in healthy "
        "plasma cells — but plasma-cell depletion is survivable, so the drug is viable. "
        "**A naive filter that just removed targets with high normal expression would have "
        "thrown BCMA away.** Which lineage carries it matters more than how much."
    )

    st.markdown("---")
    st.markdown("### Did machine learning help?")
    p = FIGURES / "story_scoreboard.png"
    if p.exists():
        st.image(str(p), width='stretch')
    st.markdown(
        "No — and reporting that is the point. Only **64 bispecific pairs** are publicly known, "
        "and **48 of them share CD3** as one arm. That is not enough data to learn from: the "
        "LightGBM model scored 0.57 on a fair test, and **0.43 — worse than guessing** — when "
        "asked about a kind of pair it had never seen. The transparent rubric, whose weights are "
        "written down and can be argued with, beat it."
    )
    st.info(
        "**The idea that makes the method work:** the two arms of a bispecific are *not "
        "interchangeable*. In a T-cell engager the tumour arm should be on the tumour, but the "
        "CD3 arm is *meant* to be on T-cells and absent from it. An earlier version of this "
        "project scored every pair the same way and came out **below chance (0.43)** — it was "
        "penalising exactly the bispecifics that work in real medicine. Splitting the scoring "
        "by design type fixed it."
    )
    st.caption(CAVEAT)

# ----------------------------------------------------------------- ranked
with tabs[1]:
    st.subheader("Ranked pairs")
    st.markdown(
        "Scores are percentiles **within design class** — a TAA×effector pair at 0.9 is "
        "not directly comparable to a TAA×TAA pair at 0.9, because the classes are "
        "judged on different evidence."
    )
    cols = [c for c in ["rank_in_class", "pair", "design_class", "rationale", "score",
                        "coverage", "is_documented", "taa_tumour_best_frac",
                        "taa_window_vital_log", "taa_normal_vital_max_log",
                        "sc_pboth_max", "sc_best_disease", "line_double_pos_frac",
                        "normal_and_max_log", "normal_double_pos_vital_n",
                        "and_gate_gain"] if c in view.columns]
    st.dataframe(
        view.sort_values("score", ascending=False)[cols].head(1000),
        width='stretch', hide_index=True,
        column_config={
            "score": st.column_config.ProgressColumn("score", min_value=0, max_value=1,
                                                     format="%.3f"),
            "coverage": st.column_config.NumberColumn("coverage", format="%.2f"),
            "sc_pboth_max": st.column_config.NumberColumn(
                "same-cell co-detect", format="%.3f",
                help="Fraction of malignant cells co-expressing both antigens"),
            "taa_tumour_best_frac": st.column_config.NumberColumn(
                "tumour positivity", format="%.2f"),
        },
    )
    st.download_button("Download this view (CSV)",
                       view[cols].to_csv(index=False).encode(),
                       "filtered_pairs.csv", "text/csv")

# ----------------------------------------------------------------- detail
with tabs[2]:
    st.subheader("Pair detail")
    if view.empty:
        st.info("No pairs match the current filters.")
    else:
        options = view.sort_values("score", ascending=False).pair.head(500).tolist()
        chosen = st.selectbox("Pair", options)
        feats, uni, normal = load_features(), load_universe(), load_normal()
        if feats.empty or uni.empty:
            st.warning("Cached feature tables not found; showing table values only.")
            st.json(view[view.pair == chosen].iloc[0].dropna().to_dict())
        else:
            from bsab import explain

            row = feats[feats.pair == chosen]
            if row.empty:
                st.warning("Pair not present in the cached feature table.")
            else:
                row = row.iloc[0]
                srow = view[view.pair == chosen].iloc[0]
                # Pull the component percentiles from scored.parquet, which has
                # the c_<rationale>__<component> columns the CSV omits.
                sc_all = load_scored()
                if not sc_all.empty:
                    hit = sc_all[sc_all.pair == chosen]
                    if not hit.empty:
                        srow = pd.concat([srow.drop(labels=[c for c in srow.index
                                                            if c in hit.columns]),
                                          hit.iloc[0]])
                normal_g = normal.reindex(list(uni.gene)).fillna(0.0)
                ev = explain.pair_evidence(row, uni, normal_g, srow)
                st.code(explain.render_text(ev), language=None)
                left, right = st.columns(2)
                with left:
                    st.markdown("**Rubric components** (percentile within class)")
                    comp = ev.get("components", {})
                    if comp:
                        st.dataframe(
                            pd.DataFrame({"component": list(comp), "percentile": list(comp.values())})
                            .sort_values("percentile", ascending=False),
                            hide_index=True, width='stretch')
                with right:
                    st.markdown("**Normal lineages expressing both** (AND-gate liability)")
                    shared = ev["safety"].get("shared_normal_cell_types") or []
                    if shared:
                        st.dataframe(pd.DataFrame(
                            shared, columns=["normal cell type", f"{ev['gene_a']} nCPM",
                                             f"{ev['gene_b']} nCPM"]),
                            hide_index=True, width='stretch')
                    else:
                        st.success("No normal cell type expresses both appreciably.")

# ----------------------------------------------------------------- indication
with tabs[3]:
    st.subheader("Per-indication target shortlists")
    by_ind = load_csv("targets_by_indication.csv")
    if by_ind.empty:
        st.info("Not generated yet.")
    else:
        st.markdown(
            "A target's value is indication-specific — CLDN18 matters in gastric cancer "
            "and essentially nowhere else. Antigens are ranked here by positivity across "
            "*that* cancer type's cell lines, discounted by normal-tissue exposure."
        )
        rec = load_csv("indication_known_recovery.csv")
        n_known = int(by_ind.in_known_bsab.sum()) if "in_known_bsab" in by_ind else 0
        if n_known == 0:
            st.warning(
                "**These shortlists recover none of the antigens already used in "
                "documented bispecifics.** That is consistent with this project's main "
                "negative finding: expression contrast penalises real targets, which "
                "typically do have substantial normal-tissue expression. Treat the "
                "shortlists as a prevalence view, not a validated ranking."
            )
        pick = st.selectbox("Cancer type", sorted(by_ind.cancer_type.unique()))
        st.dataframe(by_ind[by_ind.cancer_type == pick].drop(columns=["cancer_type"]),
                     hide_index=True, width='stretch')

# ----------------------------------------------------------------- validation
with tabs[4]:
    st.subheader("Does any of this work?")
    rv = load_csv("rubric_validation.csv")
    if not rv.empty:
        st.markdown("**Rubric vs documented pairs.** Enrichment is the ratio to the "
                    "chance rate; ROC-AUC is an upper bound because most 'negatives' "
                    "are untried rather than rejected pairs.")
        st.dataframe(rv, hide_index=True, width='stretch')
        tt = rv[rv.group.astype(str).str.contains(r"TAA\|TAA \[")]
        if len(tt) >= 2:
            st.info(
                "The two dual-antigen rationales are different hypotheses. Pathway "
                "redundancy recovers the documented pairs well; AND-gate selectivity "
                "scores **below chance** — not because the idea is wrong, but because no "
                "approved dual-TAA bispecific was designed for it, so there is no "
                "validation set. That ranking is prospective and untested."
            )
    cv = load_csv("cv_metrics.csv")
    if not cv.empty:
        st.markdown("**Supervised model under three evaluation protocols.** The random "
                    "k-fold row is included to show how much arm leakage inflates the "
                    "headline — it is not a result.")
        st.dataframe(cv, hide_index=True, width='stretch')
    pb = load_csv("prior_benchmark.csv")
    if not pb.empty:
        st.markdown("**Do hand-built expression priors beat chance at finding real "
                    "targets?** Mostly no — `enrich@k` near or below 1.0 means no "
                    "information. This is why the candidate universe stays permissive.")
        st.dataframe(pb, hide_index=True, width='stretch')
    sf = load_csv("surface_filter_audit.csv")
    if not sf.empty:
        st.markdown("**Antibody-reachability filter, audited against validated targets.**")
        st.dataframe(sf, hide_index=True, width='stretch')

# ----------------------------------------------------------------- figures
with tabs[5]:
    st.subheader("Figures")
    pngs = [p for p in sorted(FIGURES.glob("*.png")) if not p.name.startswith("story_")]
    if not pngs:
        st.info("No figures yet — they are written at the end of `run.py`.")
    for p in pngs:
        st.image(str(p), width='stretch')
        st.caption(p.stem.replace("_", " "))

# ----------------------------------------------------------------- report
with tabs[6]:
    st.subheader("Report")
    rp = RESULTS / "report.md"
    if rp.exists():
        st.markdown(rp.read_text(encoding="utf-8"))
    else:
        st.info("`results/report.md` not written yet.")
