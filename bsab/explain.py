"""Explanation generation for ranked target pairs.

Design choice: explanations are assembled from the *actual* evidence behind each
score - measured co-detection rates, named normal cell types that express both
antigens, the cancer type driving the signal - rather than written freely by a
language model. A fluent paragraph that invents a mechanism is worse than no
paragraph, because it is hard to check. So the deterministic renderer is the
default and always available.

`llm_refine()` optionally rewrites that grounded evidence into smoother prose
using OPENAI_API_KEY or ANTHROPIC_API_KEY if one is set. It is given the
structured evidence and instructed not to add facts; if no key is present the
pipeline simply uses the deterministic text, which is why nothing here depends
on network access.
"""

from __future__ import annotations

import os
import textwrap

import numpy as np
import pandas as pd

from . import config as C
from .knowledge import known_pairs_df

MECHANISM_BY_ROLE = {
    "TAA|effector": "tumour-antigen x immune-effector engager (redirects effector cells to tumour)",
    "TAA|checkpoint": "tumour-targeted checkpoint blockade (localises immune release to tumour)",
    "TAA|TAA": "dual tumour-antigen targeting (avidity gain and escape prevention)",
    "checkpoint|checkpoint": "dual checkpoint blockade",
    "checkpoint|effector": "checkpoint blockade with effector co-stimulation",
}


def _fmt_pct(x) -> str:
    return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{100 * float(x):.1f}%"


def _fmt(x, nd: int = 2) -> str:
    return "n/a" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{float(x):.{nd}f}"


def top_shared_normal_cell_types(gene_a: str, gene_b: str, normal: pd.DataFrame,
                                 n: int = 5) -> list[tuple[str, float, float]]:
    """Normal cell types with the highest min(expr_a, expr_b) - the AND-gate liability.

    These are the lineages a bispecific requiring both antigens could still
    damage, and naming them is the most actionable part of the output.
    """
    if gene_a not in normal.index or gene_b not in normal.index:
        return []
    a = normal.loc[gene_a]
    b = normal.loc[gene_b]
    both = np.minimum(a.to_numpy(), b.to_numpy())
    order = np.argsort(-both)[:n]
    return [(normal.columns[i], float(a.iloc[i]), float(b.iloc[i]))
            for i in order if both[i] > 0]


def pair_evidence(row: pd.Series, universe: pd.DataFrame, normal: pd.DataFrame,
                  scored_row: pd.Series | None = None,
                  model_contrib: pd.Series | None = None) -> dict:
    a, b = row["gene_a"], row["gene_b"]
    info = universe.set_index("gene")
    known = known_pairs_df().set_index("pair")
    pair = row["pair"]

    ev: dict[str, object] = {
        "pair": pair,
        "gene_a": a,
        "gene_b": b,
        "role_a": info.at[a, "role"] if a in info.index else "?",
        "role_b": info.at[b, "role"] if b in info.index else "?",
        "mechanism": MECHANISM_BY_ROLE.get(row.get("role_combo", ""), "unclassified pairing"),
        "is_documented": pair in known.index,
    }
    if ev["is_documented"]:
        k = known.loc[pair]
        ev["documented_as"] = f"{k['example']} ({k['status']}, {k['mechanism']})"

    ev["efficacy"] = {
        "sc_codetection_max": row.get("sc_pboth_max"),
        "sc_codetection_disease": row.get("sc_best_disease"),
        "sc_coincidence_lift": row.get("sc_lift_max"),
        "sc_datasets_covered": row.get("sc_n_datasets"),
        "cell_line_double_positive_frac": row.get("line_double_pos_frac"),
        "tumour_and_expression_log": row.get("tumour_and_max_log"),
        "best_cancer_type_a": info.at[a, "best_cancer_type"] if a in info.index else None,
        "best_cancer_type_b": info.at[b, "best_cancer_type"] if b in info.index else None,
    }
    ev["safety"] = {
        "normal_and_peak_log": row.get("normal_and_max_log"),
        "normal_double_positive_celltypes": row.get("normal_double_pos_n"),
        "vital_double_positive_celltypes": row.get("normal_double_pos_vital_n"),
        "and_gate_gain": row.get("and_gate_gain"),
        "tme_margin": row.get("sc_tme_margin"),
        "shared_normal_cell_types": top_shared_normal_cell_types(a, b, normal),
    }
    if scored_row is not None:
        ev["score"] = float(scored_row.get("score", np.nan))
        ev["coverage"] = float(scored_row.get("coverage", np.nan))
        ev["design_class"] = scored_row.get("design_class")
        ev["rationale"] = scored_row.get("rationale")
        prefix = f"c_{ev['rationale']}__" if ev.get("rationale") else "c_"
        ev["components"] = {k[len(prefix):]: round(float(v), 3)
                            for k, v in scored_row.items()
                            if k.startswith(prefix) and pd.notna(v)}
    if model_contrib is not None:
        top = model_contrib.drop(labels=["_bias"], errors="ignore")
        top = top.reindex(top.abs().sort_values(ascending=False).index)[:6]
        ev["model_drivers"] = {k: round(float(v), 4) for k, v in top.items()}
    return ev


def render_text(ev: dict) -> str:
    a, b = ev["gene_a"], ev["gene_b"]
    eff, saf = ev["efficacy"], ev["safety"]
    lines: list[str] = []

    head = f"{a} x {b} - {ev['mechanism']}"
    lines.append(head)
    lines.append("-" * len(head))

    if ev["is_documented"]:
        lines.append(f"STATUS: already documented - {ev['documented_as']}. "
                     f"Recovered by the ranking rather than proposed by it.")
    else:
        lines.append("STATUS: not among the documented bispecific pairs in this project's "
                     "reference list - i.e. a genuine suggestion, with all the caveats below.")

    if "score" in ev:
        rat = ev.get("rationale")
        rat_txt = f", scored as {str(rat).replace('_', ' ')}" if rat else ""
        lines.append(f"Within-class score {_fmt(ev['score'], 3)} "
                     f"(percentile among {ev.get('design_class', 'its class')} pairs"
                     f"{rat_txt}; evidence coverage {_fmt_pct(ev.get('coverage'))} "
                     f"of rubric weight).")
        if ev.get("rationale") == "and_gate_selectivity":
            lines.append("  NOTE: the AND-gate rationale has no historical validation set - "
                         "no approved dual-antigen bispecific was designed for it - so this "
                         "ranking is a prospective hypothesis, not a validated prediction.")

    # --- efficacy ---
    cov = eff.get("sc_datasets_covered")
    if cov is not None and not (isinstance(cov, float) and np.isnan(cov)) and cov > 0:
        lines.append(
            f"EFFICACY: in single tumour cells, {_fmt_pct(eff['sc_codetection_max'])} of malignant "
            f"cells co-express both antigens (strongest in {eff['sc_codetection_disease']}; "
            f"measured across {int(cov)} dataset(s)). Co-incidence lift "
            f"{_fmt(eff['sc_coincidence_lift'])}x versus independent expression"
            + (" - the two antigens tend to occupy the same cells."
               if (eff['sc_coincidence_lift'] or 0) > 1.2 else
               " - little evidence of same-cell preference.")
        )
    else:
        lines.append("EFFICACY: no single-cell co-detection available (one or both genes absent "
                     "from the tumour datasets used), so same-cell co-expression is unverified.")
    lines.append(
        f"  Across 1,206 cancer cell lines, {_fmt_pct(eff['cell_line_double_positive_frac'])} "
        f"carry both at nTPM>=10. Individually the arms peak in "
        f"{eff['best_cancer_type_a']} ({a}) and {eff['best_cancer_type_b']} ({b})."
    )

    # --- safety ---
    shared = saf.get("shared_normal_cell_types") or []
    vital = saf.get("vital_double_positive_celltypes")
    lines.append(
        f"SAFETY: the AND-gate peak in normal tissue is {_fmt(saf['normal_and_peak_log'])} "
        f"(log1p nCPM); {_fmt(saf['normal_double_positive_celltypes'], 0)} of 154 normal cell types "
        f"express both above threshold, of which {_fmt(vital, 0)} are dose-limiting lineages. "
        f"Pairing lowers peak normal exposure by {_fmt(saf['and_gate_gain'])} log units versus the "
        f"worse single arm."
    )
    if shared:
        named = "; ".join(f"{ct} ({a} {va:.0f}, {b} {vb:.0f} nCPM)" for ct, va, vb in shared[:4])
        lines.append(f"  Normal lineages at risk (highest shared expression): {named}.")
    else:
        lines.append("  No normal cell type shows appreciable expression of both antigens.")
    tme = saf.get("tme_margin")
    if tme is not None and not (isinstance(tme, float) and np.isnan(tme)):
        lines.append(
            f"  Within the tumour itself, double-positive rate is {_fmt_pct(tme)} higher in "
            f"malignant than in non-malignant cells"
            + (" - a favourable in-situ margin." if tme > 0 else
               " - i.e. NOT selective for malignant cells in these samples.")
        )

    if "model_drivers" in ev:
        drivers = ", ".join(f"{k} {v:+.3f}" for k, v in ev["model_drivers"].items())
        lines.append(f"MODEL: top additive contributions - {drivers}.")

    lines.append(
        "CAVEATS: expression is mRNA, not surface protein; scRNA-seq dropout makes "
        "co-detection an under-estimate; cell lines lack microenvironment and often lose "
        "lineage markers. This is target triage to prioritise experiments, not evidence of efficacy."
    )
    return "\n".join(textwrap.fill(l, 100, subsequent_indent="  ") if not set(l) <= {"-"} else l
                     for l in lines)


# --------------------------------------------------------------------------
# Optional LLM rewrite
# --------------------------------------------------------------------------
LLM_SYSTEM = (
    "You rewrite structured evidence about an antibody target pair into a short, sober "
    "paragraph for a drug-discovery audience. Rules: use ONLY the numbers and facts given; "
    "never introduce a mechanism, citation, trial, or claim that is not in the input; keep "
    "every caveat; do not express enthusiasm. If the evidence is weak, say so plainly."
)


def llm_available() -> str | None:
    if os.environ.get("ANTHROPIC_API_KEY"):
        return "anthropic"
    if os.environ.get("OPENAI_API_KEY"):
        return "openai"
    return None


def llm_refine(ev: dict, deterministic: str, model: str | None = None) -> tuple[str, str]:
    """Rewrite the grounded text with an LLM if a key is configured.

    Returns (text, provenance). Any failure falls back to the deterministic
    text, so the pipeline never depends on this succeeding.
    """
    provider = llm_available()
    if provider is None:
        return deterministic, "deterministic (no LLM key configured)"
    import json

    payload = json.dumps(ev, default=str, indent=2)
    prompt = (f"Structured evidence:\n{payload}\n\n"
              f"Deterministic summary to rewrite:\n{deterministic}")
    try:
        if provider == "anthropic":
            import anthropic

            client = anthropic.Anthropic()
            resp = client.messages.create(
                model=model or "claude-opus-5-5",
                max_tokens=700,
                system=LLM_SYSTEM,
                messages=[{"role": "user", "content": prompt}],
            )
            return resp.content[0].text.strip(), f"anthropic:{model or 'claude-opus-5-5'}"
        from openai import OpenAI

        client = OpenAI()
        resp = client.chat.completions.create(
            model=model or "gpt-4o",
            messages=[{"role": "system", "content": LLM_SYSTEM},
                      {"role": "user", "content": prompt}],
            max_tokens=700,
        )
        return resp.choices[0].message.content.strip(), f"openai:{model or 'gpt-4o'}"
    except Exception as exc:  # noqa: BLE001 - explanation is non-critical
        return deterministic, f"deterministic (LLM call failed: {type(exc).__name__})"


def explain_top(features: pd.DataFrame, scored: pd.DataFrame, universe: pd.DataFrame,
                normal: pd.DataFrame, n: int = 25, use_llm: bool = False,
                model_contrib: pd.DataFrame | None = None) -> list[dict]:
    carry = ["pair", "score", "coverage", "design_class", "rationale"]
    carry = [c for c in carry if c in scored.columns]
    merged = features.merge(scored[carry + [c for c in scored.columns if c.startswith("c_")]],
                            on="pair", how="left")
    # Take the best from each design class in turn, so a class with many pairs
    # cannot crowd the others out of the explanation set.
    merged = merged.sort_values("score", ascending=False)
    if "design_class" in merged.columns:
        per_class = max(1, n // max(merged.design_class.nunique(), 1))
        top = (merged.groupby("design_class", group_keys=False)
               .head(per_class).sort_values("score", ascending=False).head(n))
        if len(top) < n:
            top = pd.concat([top, merged[~merged.pair.isin(top.pair)].head(n - len(top))])
    else:
        top = merged.head(n)
    out = []
    for _, row in top.iterrows():
        keep = [c for c in row.index if c.startswith("c_")]
        keep += [c for c in ("score", "coverage", "design_class", "rationale") if c in row.index]
        srow = row[keep]
        contrib = None
        if model_contrib is not None and row["pair"] in model_contrib.index:
            contrib = model_contrib.loc[row["pair"]]
        ev = pair_evidence(row, universe, normal, srow, contrib)
        text = render_text(ev)
        if use_llm:
            text, prov = llm_refine(ev, text)
        else:
            prov = "deterministic"
        out.append({"pair": row["pair"], "score": row["score"],
                    "evidence": ev, "text": text, "provenance": prov})
    return out
