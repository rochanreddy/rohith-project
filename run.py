#!/usr/bin/env python
"""Run the full bispecific target-pair pipeline.

    python run.py                 # full run (cached stages are reused)
    python run.py --force         # recompute every stage
    python run.py --llm           # rewrite explanations with an LLM if a key is set
    python run.py --top 40        # number of pairs to explain in detail
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from bsab import config as C
from bsab import explain, features_pair, model, pipeline, report, score



def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="recompute all cached stages")
    ap.add_argument("--llm", action="store_true", help="LLM-rewrite explanations if a key is set")
    ap.add_argument("--top", type=int, default=25, help="pairs to explain in detail")
    ap.add_argument("--skip-model", action="store_true", help="rubric only, no supervised model")
    args = ap.parse_args()

    ctx = pipeline.run_all(force=args.force)
    universe, pairs, feats, scored = ctx["universe"], ctx["pairs"], ctx["features"], ctx["scored"]
    labels = features_pair.label_pairs(feats[["pair"]])

    ctx["rubric_components"] = score.component_table()
    ctx["prior_benchmark"] = (pd.read_csv(C.RESULTS / "prior_benchmark.csv")
                              if (C.RESULTS / "prior_benchmark.csv").exists() else pd.DataFrame())
    ctx["sc_summary"] = (pd.read_csv(C.RESULTS / "scrna_datasets.csv")
                         if (C.RESULTS / "scrna_datasets.csv").exists() else pd.DataFrame())
    from bsab import targets as T

    ctx["surface_audit"] = T.audit_surface_filter(verbose=False)

    # ---------------- ranked output ----------------
    pipeline._stamp("writing ranked pair tables")
    ranked = (feats.merge(scored.drop(columns=["pair"]).assign(pair=scored["pair"].to_numpy()),
                          on="pair", how="left")
              # bool, not int: `~` on an int column is a bitwise NOT, which
              # pandas then reads as a column selector rather than a mask.
              .assign(is_documented=labels.to_numpy().astype(bool))
              .sort_values(["design_class", "score"], ascending=[True, False]))
    ranked["rank_in_class"] = ranked.groupby("design_class").cumcount() + 1
    ranked = ranked.sort_values("score", ascending=False)
    ranked.insert(0, "rank_overall", range(1, len(ranked) + 1))
    # The full ~200-column matrix is 240 MB as CSV, so it goes to parquet and the
    # CSV carries the columns a reader actually scans.
    ranked.to_parquet(C.RESULTS / "ranked_pairs_full.parquet", index=False)

    headline = ["rank_overall", "rank_in_class", "pair", "design_class", "score", "coverage",
                "low_confidence", "is_documented", "taa_surface_confidence",
                "surface_confidence__min", "taa_tumour_best_frac",
                "taa_window_vital_log", "taa_normal_vital_max_log", "sc_pboth_max",
                "sc_best_disease", "sc_lift_max", "line_double_pos_frac",
                "normal_and_max_log", "normal_double_pos_vital_n", "and_gate_gain"]
    headline = [c for c in headline if c in ranked.columns]
    score_cols = [c for c in ranked.columns if c.startswith(("score__", "coverage__"))]
    ranked[headline + score_cols].to_csv(C.RESULTS / "ranked_pairs.csv", index=False)
    ranked[ranked.is_documented][headline].to_csv(
        C.RESULTS / "documented_pairs_ranked.csv", index=False)

    # Per-class candidate tables: scores are within-class percentiles, so the
    # useful output is the top of each class, not one blended list.
    #
    # TAA x anchor scores depend only on the tumour arm, so every anchor partner
    # for a given antigen ties exactly. Listing all of them would fill the table
    # with near-duplicates, so the headline keeps the best partner per tumour
    # antigen and `n_anchor_partners_tied` records how many tied with it. The
    # full ranking in ranked_pairs.csv keeps every row.
    cand = ranked[~ranked.is_documented].copy()
    cand["taa_arm"] = np.where(
        cand.gene_a.isin(set(universe.gene[universe.role == "TAA"])), cand.gene_a, cand.gene_b)
    novel_parts = []
    for cls, grp in cand.groupby("design_class"):
        if cls in ("TAA|effector", "TAA|checkpoint"):
            tie = grp.groupby(["taa_arm", "score"]).pair.transform("size")
            grp = grp.assign(n_anchor_partners_tied=tie)
            grp = grp.sort_values("rank_in_class").drop_duplicates("taa_arm")
        novel_parts.append(grp.nsmallest(40, "rank_in_class"))
    extra = ["n_anchor_partners_tied", "surface_confidence__min"]
    out_cols = headline + [c for c in extra if any(c in p.columns for p in novel_parts)]
    novel = (pd.concat(novel_parts)
             .reindex(columns=[c for c in out_cols if c in pd.concat(novel_parts).columns])
             .sort_values(["design_class", "rank_in_class"]))
    novel.to_csv(C.RESULTS / "novel_candidates.csv", index=False)
    ctx["novel_top"] = novel
    ctx["ranked"] = ranked

    # ---------------- per-indication shortlists ----------------
    pipeline._stamp("building per-indication target shortlists")
    from bsab import indication

    by_ind = indication.rank_targets_by_indication(universe)
    by_ind.to_csv(C.RESULTS / "targets_by_indication.csv", index=False)
    ind_recall = indication.recovered_known_targets(by_ind)
    ind_recall.to_csv(C.RESULTS / "indication_known_recovery.csv", index=False)
    ctx["by_indication"] = by_ind
    ctx["indication_recovery"] = ind_recall
    print(f"    {by_ind.cancer_type.nunique()} indications, "
          f"{int(by_ind.in_known_bsab.sum())} documented targets recovered in shortlists")

    # ---------------- supervised model (secondary) ----------------
    cv_rows: list[dict] = []
    importance = pd.DataFrame()
    contrib = None
    if not args.skip_model:
        cv_path = C.RESULTS / "cv_metrics.csv"
        if cv_path.exists() and not args.force:
            cv_rows = pd.read_csv(cv_path).to_dict("records")
            pipeline._stamp(f"reusing cached cross-validation ({cv_path.name})")
        else:
            pipeline._stamp("supervised model: arm-grouped cross-validation")
            grouped = model.cross_validate(feats, pairs, labels, exclude_cd3e=False)
            cv_rows.append(grouped["metrics"])
            print("   ", grouped["metrics"])

            pipeline._stamp("supervised model: CD3E removed from training")
            no_cd3e = model.cross_validate(feats, pairs, labels, exclude_cd3e=True)
            cv_rows.append(no_cd3e["metrics"])
            print("   ", no_cd3e["metrics"])

            pipeline._stamp("supervised model: naive random k-fold (for contrast)")
            leaky = model.random_cv_for_contrast(feats, labels)
            cv_rows.append(leaky["metrics"])
            print("   ", leaky["metrics"])
            pd.DataFrame(cv_rows).to_csv(cv_path, index=False)

        pipeline._stamp("fitting final model for feature importance")
        clf, importance, cols = model.fit_full(feats, labels)
        importance.to_csv(C.RESULTS / "feature_importance.csv", index=False)



        # SHAP contributions only for pairs that might be explained - computing
        # them for all 232k rows costs minutes and nothing reads the rest.
        wanted = set(ranked.groupby("design_class").head(100).pair)
        rows = np.flatnonzero(feats.pair.isin(wanted).to_numpy())
        contrib = model.contributions(clf, feats, cols, rows=rows)
        contrib.index = feats.pair.to_numpy()[rows]
    ctx["cv_rows"] = cv_rows
    ctx["importance"] = importance

    # ---------------- explanations ----------------
    pipeline._stamp(f"generating explanations for top {args.top} pairs")
    # Raw nCPM (not log) so explanations can quote real expression values.
    from bsab import data_hpa as D

    normal_raw = D.load_normal_sc().reindex(list(universe.gene)).fillna(0.0)
    expls = explain.explain_top(feats, scored, universe, normal_raw,
                                n=args.top, use_llm=args.llm, model_contrib=contrib)
    ctx["explanations"] = expls
    with (C.RESULTS / "explanations.md").open("w", encoding="utf-8") as fh:
        fh.write("# Explanations for top-ranked bispecific target pairs\n\n")
        fh.write(f"Provenance: {expls[0]['provenance'] if expls else 'n/a'}\n\n")
        for e in expls:
            fh.write("```\n" + e["text"] + "\n```\n\n")
    if args.llm and explain.llm_available() is None:
        print("    note: --llm requested but no API key configured; used deterministic text")

    # ---------------- figures ----------------
    pipeline._stamp("rendering figures")
    report.fig_prior_benchmark(ctx["prior_benchmark"], C.FIGURES / "prior_benchmark.png")
    report.fig_cv_comparison(cv_rows, C.FIGURES / "cv_comparison.png")
    if not importance.empty:
        report.fig_feature_importance(importance, C.FIGURES / "feature_importance.png")
    report.fig_safety_efficacy(feats, scored, C.FIGURES / "safety_efficacy.png")
    report.fig_sc_heatmap(feats, scored, C.FIGURES / "sc_codetection.png")

    pipeline._stamp("writing report")
    report.write_report(ctx)

    print("\nDone. Key outputs:")
    for p in ["report.md", "ranked_pairs.csv", "novel_candidates.csv", "explanations.md",
              "cv_metrics.csv", "prior_benchmark.csv", "rubric_validation.csv"]:
        f = C.RESULTS / p
        if f.exists():
            print(f"  results/{p}  ({f.stat().st_size / 1024:.0f} KB)")
    print(f"  results/figures/  ({len(list(C.FIGURES.glob('*.png')))} figures)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
