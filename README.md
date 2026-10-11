# Bispecific antibody target-pair prediction

Ranks two-target combinations for bispecific-antibody development in oncology,
using public single-cell and tumour expression data, biological feature
engineering, pairwise learning, and grounded explanation generation.

The central question: given a few thousand antibody-accessible proteins, which
*pair* is worth laboratory work? Efficacy needs both antigens on the **same**
tumour cell. Safety is set by whichever normal tissue also carries both — and
for a bispecific requiring dual engagement, exposure in a tissue scales with the
**lesser** of the two antigens, not either alone. That asymmetry is the point of
pairing, and it is what this pipeline measures.

## What it produces

| Output | Content |
|---|---|
| `results/report.md` | Full written report: method, findings, limitations |
| `results/ranked_pairs.csv` | All 187,796 pairs ranked (headline columns) |
| `results/ranked_pairs_full.parquet` | Same, with all ~200 features |
| `results/novel_candidates.csv` | Top-ranked pairs not in the reference list |
| `results/explanations.md` | Evidence-grounded rationale per top pair |
| `results/cv_metrics.csv` | Model evaluation under three protocols |
| `results/prior_benchmark.csv` | Do expression heuristics beat chance? |
| `results/surface_filter_audit.csv` | Reachability filter vs validated targets |
| `results/figures/` | Five figures (CSV table beside each) |

## Quick start

```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements.txt     # Linux/macOS: .venv/bin/pip

python scripts/download_data.py                   # ~1.4 GB, public sources
python run.py                                     # full pipeline
.venv/Scripts/streamlit run app.py                # browse results -> localhost:8501
```

The Streamlit app reads whatever `run.py` produced: filter and sort the ranked
pairs, open any pair to see its evidence-grounded explanation and the named
normal lineages that express both antigens, browse per-indication shortlists,
and read the validation tables and report.

Useful flags: `--force` recompute all cached stages, `--top N` pairs to explain,
`--skip-model` rubric only, `--llm` rewrite explanations with a language model
(needs `OPENAI_API_KEY` or `ANTHROPIC_API_KEY`; falls back silently if absent).

Stages cache to `data/processed/`, so reruns resume. A cold run is a few minutes
apart from the single-cell stage (~6 min) and the 1 GB cell-line pivot (~2 min).

## Data

| Source | Content | Role |
|---|---|---|
| HPA single-cell atlas | 20,151 genes × 154 normal cell types | normal-tissue safety |
| HPA cell lines | 20,151 genes × 1,206 lines, 25 cancer types | tumour expression |
| HPA annotations | protein class, localisation, secretome | antibody accessibility |
| CELLxGENE | 5 tumour scRNA-seq datasets, 47,061 malignant cells | same-cell co-expression |
| `bsab/knowledge.py` | 64 documented bispecific pairs | positive labels |

## Method

1. **Universe** — antibody-reachable genes with some tumour expression: 3,852.
   Reachability is *graded* (high / medium / none), not a single flag — see below.
2. **Pair space** — bounded *structurally*, not by a quality prior: every
   tumour antigen × every immune anchor, anchor × anchor, and TAA × TAA among
   the 400 most prevalent antigens. 187,796 pairs, containing all 64 known ones.
3. **Features** — ~200 per pair: symmetric aggregates (min/max/mean/|diff|) of 37
   single-arm features so arm order cannot be encoded, features resolved to the
   *tumour arm* specifically, plus pairwise quantities (AND-gate peaks in normal
   and tumour tissue, double-positive counts, profile correlations,
   functional-term overlap, single-cell co-detection).
4. **Ranking** — declared percentile rubrics with published weights, **split by
   design class** because the arms of a bispecific are not interchangeable.
   Components lacking data are dropped and remaining weights renormalised, so
   unmeasured pairs are not penalised; `coverage` records how much was available.
5. **Supervised model** — LightGBM, reported as a *secondary* signal with
   arm-grouped cross-validation and a CD3E-excluded generalisation test.
6. **Explanations** — assembled deterministically from the evidence behind each
   score, including the named normal cell types expressing both antigens.

## Three findings that shaped the design

**Expression heuristics barely beat chance.** Four hand-built single-target
priors were benchmarked against the documented-bispecific list. The best reaches
~2× enrichment; ranking by tumour prevalence, or by *low* normal expression,
does **worse than random**. Validated targets like EGFR, EPCAM, MS4A1 and MET
have substantial normal expression — what makes them usable is whether the
normal lineage is expendable and whether toxicity is clinically manageable, which
this data cannot see. A pipeline triaging on expression contrast would have
discarded most drugs that work. Hence the permissive universe.

**A symmetric rubric is wrong, and it scored below chance.** The first version
rewarded same-cell co-expression for every pair and reached ROC-AUC **0.43** —
worse than random. The reason is biological: in a T-cell engager the effector arm
(CD3E) is *supposed* to be on T cells and absent from tumour, so same-cell
co-detection is near zero by design. The metric was penalising the dominant class
of real bispecifics. Scoring is now per design class, and within dual-TAA pairs
the two real rationales separate sharply — pathway redundancy reaches ROC-AUC
0.83 while AND-gate selectivity sits at 0.30, because no approved dual-TAA
bispecific was *designed* for AND-gating. That ranking is prospective and
untested, not validated.

Relatedly: HPA's "Predicted membrane proteins" is a topology call that includes
every organelle membrane. Using it directly filled the top of the rankings with
mitochondrial ATP-synthase subunits, ER enzymes (VKORC1 — the warfarin target)
and lysosomal IFI30. Reachability is now graded, and the filter is audited
against the validated targets: 46 high, 4 medium, 6 rejected of 56.

**The labels are too few and too lopsided for ML to be the headline.** 64
documented pairs, **48 sharing CD3E** as one arm, and unlabelled pairs are
untried rather than rejected (positive-unlabelled). A model trained this way
largely learns "one arm should be CD3E". Worse, random k-fold puts `CD19|CD3E`
in training and `CD22|CD3E` in test, inflating the headline substantially — the
pipeline reports that leaky number alongside the honest one precisely to show the
gap. So the transparent rubric is primary and the model is corroborating.

## Limitations

mRNA is not surface protein. scRNA-seq dropout makes co-detection a systematic
under-estimate. Cell lines lack microenvironment and lose lineage markers
(CLDN18, SSTR2, ANGPT2, DLL4 sit near the floor despite being real targets).
Only five tumour datasets. The rubric weights are asserted, not fitted. No
immunology beyond expression — effector fitness, cytokine-release risk, antigen
density thresholds, resistance and manufacturability all decide programmes and
are out of scope. The universe cut-off was set knowing the reference list, so
in-universe recall is optimistic; trust the arm-grouped evaluations.

**This is target triage to prioritise experiments. It is not evidence of
efficacy or safety, and nothing here should inform a clinical decision.**

## Layout

```
bsab/
  knowledge.py      64 curated bispecific pairs (the labels)
  config.py         paths, thresholds, dose-limiting lineages
  data_hpa.py       HPA loaders with parquet caching
  targets.py        universe, prior benchmark, pair space
  features_single.py  per-arm features
  features_sc.py    single-cell co-detection
  pairops.py        chunked pairwise reductions
  features_pair.py  pair feature assembly
  score.py          the transparent rubric
  model.py          LightGBM + leakage-aware evaluation
  explain.py        grounded explanations (+ optional LLM)
  report.py         figures and report
  pipeline.py       cached stages
run.py              CLI entry point
scripts/download_data.py
```
