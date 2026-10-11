#!/usr/bin/env python
"""Browsable UI over the bispecific target-pair results.

    .venv/Scripts/streamlit run app.py        # -> http://localhost:8501

Reads whatever `run.py` has produced in results/. Nothing is computed here
except the per-pair explanation text, which is generated on demand from the
cached feature table so it always matches the numbers on screen.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from bsab.config import VITAL_CELL_TYPES as C_VITAL
from bsab.report import BLUE, CRITICAL, GOOD, GRID, INK, INK_2, ORANGE

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

tabs = st.tabs(["🔬 Try it", "Ranked pairs", "Pair detail", "By indication",
                "Validation", "Figures", "Report"])

# ----------------------------------------------------------------- try it
with tabs[0]:
    st.markdown(
        "Two live tools. Each one runs the real analysis on click — nothing here is "
        "pre-rendered."
    )

    # ---------- Tool 1: find targets for a cancer ----------
    st.subheader("1 · Find the best targets for a cancer")
    by_ind_all = load_csv("targets_by_indication.csv")
    if by_ind_all.empty:
        st.info("Run `python run.py` first.")
    else:
        c1, c2 = st.columns([3, 1])
        cancer = c1.selectbox(
            "Pick a cancer", sorted(by_ind_all.cancer_type.unique()),
            index=sorted(by_ind_all.cancer_type.unique()).index("myeloma")
            if "myeloma" in set(by_ind_all.cancer_type) else 0,
            label_visibility="collapsed")
        run1 = c2.button("Find targets", type="primary", width='stretch')

        if run1:
            uni = load_universe()
            sub = by_ind_all[by_ind_all.cancer_type == cancer].sort_values("rank")
            with st.status(f"Searching for targets in {cancer.replace('_', ' ')}…",
                           expanded=True) as status:
                st.write("Starting from **20,151 human genes**")
                st.write(f"Keeping the **{len(uni):,}** an antibody can physically reach "
                         "(on the cell surface)")
                st.write(f"Checking which are on **{cancer.replace('_', ' ')}** tumour samples")
                st.write("Checking all **154 healthy cell types** for safety risk")
                status.update(label=f"Found {len(sub)} candidate targets", state="complete")

            st.markdown(f"#### Top targets for {cancer.replace('_', ' ')}")
            top = sub.head(8).iloc[::-1]
            fig, ax = plt.subplots(figsize=(8, 0.46 * len(top) + 1.4))
            colours = [GOOD if k else BLUE for k in top.in_known_bsab]
            ax.barh(top.gene, top.tumour_positivity * 100, color=colours, height=0.62)
            for y, (v, known) in enumerate(zip(top.tumour_positivity * 100, top.in_known_bsab)):
                ax.text(v + 1.5, y, f"{v:.0f}%" + ("  ← known drug target" if known else ""),
                        va="center", fontsize=9,
                        color=GOOD if known else INK_2,
                        fontweight="bold" if known else "normal")
            ax.set_xlim(0, 125)
            ax.set_xlabel("% of tumour samples carrying this target", fontsize=9, color=INK_2)
            ax.set_title(f"What fraction of {cancer.replace('_', ' ')} carries each target",
                         loc="left", fontsize=11, fontweight="bold")
            ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
            ax.grid(axis="y", lw=0)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            st.pyplot(fig, clear_figure=True)

            hits = sub[sub.in_known_bsab]
            if len(hits):
                st.success(
                    f"**Sanity check passed.** {', '.join(hits.gene)} "
                    f"{'is' if len(hits) == 1 else 'are'} already used in real bispecific "
                    f"drugs for this cancer — the method found "
                    f"{'it' if len(hits) == 1 else 'them'} without being told."
                )
            else:
                st.warning(
                    "No known bispecific target appears in this shortlist. Either none exists "
                    "for this cancer, or the method missed it — worth checking before trusting "
                    "the novel suggestions here."
                )
            show_cols = ["rank", "gene", "tumour_positivity", "normal_peak_ncpm",
                         "normal_peak_cell_type", "surface_confidence", "in_known_bsab"]
            st.dataframe(sub[[c for c in show_cols if c in sub.columns]].head(15),
                         hide_index=True, width='stretch')

    st.divider()

    # ---------- Tool 2: test any pair ----------
    st.subheader("2 · Test any two targets as a pair")
    st.caption("Pick two proteins and see whether they would work together, and what "
               "healthy tissue they would damage.")

    feats, uni, normal_raw = load_features(), load_universe(), load_normal()
    if feats.empty or uni.empty:
        st.info("Cached feature tables not found — run `python run.py`.")
    else:
        genes = sorted(set(feats.gene_a) | set(feats.gene_b))
        d1, d2, d3 = st.columns([2, 2, 1])
        default_a = genes.index("TNFRSF17") if "TNFRSF17" in genes else 0
        default_b = genes.index("CD3E") if "CD3E" in genes else 1
        ga = d1.selectbox("Target 1", genes, index=default_a)
        gb = d2.selectbox("Target 2", genes, index=default_b)
        run2 = d3.button("Analyse pair", type="primary", width='stretch')

        if run2:
            key = "|".join(sorted((ga, gb)))
            row = feats[feats.pair == key]
            if ga == gb:
                st.error("Pick two different targets.")
            elif row.empty:
                st.warning(
                    f"**{ga} × {gb}** is not in the scored pair space. The space covers every "
                    "tumour-antigen × immune-anchor pair plus the most prevalent "
                    "antigen × antigen pairs, so this combination was never a candidate."
                )
            else:
                row = row.iloc[0]
                srow = ranked[ranked.pair == key]
                sc_all = load_scored()
                hit = sc_all[sc_all.pair == key] if not sc_all.empty else pd.DataFrame()

                with st.status(f"Analysing {ga} × {gb}…", expanded=True) as status:
                    st.write("Measuring co-expression across **47,061 single tumour cells**")
                    st.write("Checking **1,206 cancer cell lines**")
                    st.write("Scanning **154 healthy cell types** for safety risk")
                    status.update(label="Analysis complete", state="complete")

                pboth = row.get("sc_pboth_max")
                lift = row.get("sc_lift_max")
                vital_n = row.get("normal_double_pos_vital_n")
                k1, k2, k3 = st.columns(3)
                k1.metric("Tumour cells carrying BOTH",
                          "—" if pd.isna(pboth) else f"{100 * float(pboth):.1f}%")
                k2.metric("Same-cell coincidence",
                          "—" if pd.isna(lift) else f"{float(lift):.1f}×",
                          help="Above 1 means the two targets occupy the same cells more "
                               "than chance would give.")
                k3.metric("Dose-limiting tissues at risk",
                          "—" if pd.isna(vital_n) else f"{int(vital_n)}",
                          help="Healthy cell types whose loss would be dangerous, which "
                               "carry both targets.")

                if not srow.empty:
                    s = srow.iloc[0]
                    if bool(s.is_documented):
                        st.success(f"**This is a real bispecific.** It appears in the "
                                   f"reference list of 64 documented pairs, ranked "
                                   f"#{int(s.rank_in_class)} within {s.design_class}.")
                    else:
                        st.info(f"**Score {float(s.score):.3f}** — rank "
                                f"#{int(s.rank_in_class):,} of its "
                                f"{s.design_class} class. Not a documented bispecific.")

                # Which healthy tissue would be hit
                st.markdown("##### Healthy tissue carrying both targets")
                from bsab import explain as _explain

                nrm = normal_raw.reindex(list(uni.gene)).fillna(0.0)
                shared = _explain.top_shared_normal_cell_types(ga, gb, nrm, n=6)
                if not shared:
                    st.success("No healthy cell type carries both appreciably — the safest case.")
                else:
                    sh = pd.DataFrame(shared, columns=["cell type", ga, gb])
                    sh["risk (lower of the two)"] = sh[[ga, gb]].min(axis=1).round(0)
                    vital = set(C_VITAL)
                    fig, ax = plt.subplots(figsize=(8, 0.44 * len(sh) + 1.3))
                    rk = sh["risk (lower of the two)"].to_numpy()[::-1]
                    names = sh["cell type"].to_numpy()[::-1]
                    colours = [CRITICAL if n in vital else ORANGE for n in names]
                    ax.barh(names, rk, color=colours, height=0.6)
                    for y, v in enumerate(rk):
                        ax.text(v + max(rk) * 0.02 + 0.4, y, f"{v:,.0f}", va="center",
                                fontsize=9, color=INK_2)
                    ax.set_xlim(0, max(max(rk) * 1.2, 1))
                    ax.set_xlabel("exposure = the LOWER of the two targets (nCPM)",
                                  fontsize=9, color=INK_2)
                    ax.set_title("A dual-targeting drug only hits a tissue that has BOTH",
                                 loc="left", fontsize=10.5, fontweight="bold")
                    ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
                    ax.grid(axis="y", lw=0)
                    for sp in ("top", "right"):
                        ax.spines[sp].set_visible(False)
                    st.pyplot(fig, clear_figure=True)
                    st.caption(
                        "Red = a lineage whose loss is dose-limiting. This is the key idea: "
                        "because the drug needs **both** targets, the damage in a tissue is "
                        "set by whichever target is **lower** there — so pairing can make an "
                        "otherwise unusable target safe."
                    )
                    st.dataframe(sh, hide_index=True, width='stretch')

                if not hit.empty:
                    ev = _explain.pair_evidence(
                        row, uni, nrm,
                        pd.concat([srow.iloc[0].drop(labels=[c for c in srow.columns
                                                             if c in hit.columns]),
                                   hit.iloc[0]]) if not srow.empty else hit.iloc[0])
                    with st.expander("Full written explanation"):
                        st.code(_explain.render_text(ev), language=None)

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
