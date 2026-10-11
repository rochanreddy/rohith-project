"""Definition of the candidate target universe and the pair space.

A bispecific arm has to be reachable by an antibody, so the universe starts
from HPA's predicted membrane / secreted proteins plus proteins localised to
the plasma membrane (~8,000 genes), and keeps those with some tumour expression.

We deliberately do **not** shrink that universe with a hand-built "target
quality" prior. `benchmark_priors()` measures several such priors against the
documented-bispecific list and finds all of them close to chance, so using one
to pick the universe would discard validated targets while pretending to focus
the search. Instead the pair space is bounded structurally:

  * every tumour antigen x every immune anchor (effector receptor / checkpoint)
  * anchor x anchor, excluding effector x effector
  * tumour antigen x tumour antigen, restricted to the `N_TOP_TAA_FOR_TAA_PAIRS`
    most prevalent antigens by tumour-line positivity

That mirrors how real bispecifics are actually built and needs only one
transparent cut-off.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C
from . import data_hpa as D
from .knowledge import CHECKPOINT_GENES, EFFECTOR_GENES, known_pairs_df

MEMBRANE_CLASS = "Predicted membrane proteins"
SECRETED_CLASS = "Predicted secreted proteins"


INTRACELLULAR_CLASS = "Predicted intracellular proteins"

# Compartments that place a protein out of an antibody's reach. HPA's
# "Predicted membrane proteins" class is a *topology* call - it fires for every
# organelle membrane - so without this exclusion the rankings fill up with
# mitochondrial ATP-synthase subunits, ER enzymes and lysosomal proteins that
# no antibody can bind in a living cell.
ORGANELLE_TERMS = (
    "Mitochondria|Endoplasmic reticulum|Golgi|Nucleo|Nuclear|Cytosol|Lysosome|"
    "Peroxisome|Vesicles|Cytoskeleton|Centrosome|Actin|Microtubule|Mitotic|"
    "Intermediate filaments|Rods|Aggresome|Nucleoli"
)


def _accessibility(annot: pd.DataFrame) -> pd.DataFrame:
    """Grade how reachable each protein is by an antibody.

    `surface_confidence` is "high" when there is positive evidence of surface or
    extracellular residence (plasma-membrane localisation, secretion, CD marker,
    GPCR), "medium" for a predicted membrane protein with no contradicting
    organelle evidence, and "none" otherwise.

    Graded rather than hard-filtered on purpose: requiring plasma-membrane
    evidence alone would drop ten validated bispecific targets (CLDN18, FOLH1,
    ROR1, STEAP1, TPBG, PMEL, GUCY2C, DLL4, CLDN6, TIGIT), which simply lack
    reliable immunofluorescence data in HPA. `audit_surface_filter()` reports
    that false-negative rate instead of hiding it.
    """
    pclass = annot["Protein class"].map(D.split_terms)
    is_membrane = pclass.map(lambda s: MEMBRANE_CLASS in s).fillna(False)
    is_secreted = pclass.map(lambda s: SECRETED_CLASS in s).fillna(False)
    is_intracellular = pclass.map(lambda s: INTRACELLULAR_CLASS in s).fillna(False)
    is_cd_marker = pclass.map(lambda s: "CD markers" in s).fillna(False)
    is_gpcr = pclass.map(lambda s: "G-protein coupled receptors" in s).fillna(False)
    is_tf = pclass.map(lambda s: "Transcription factors" in s).fillna(False)
    is_ribosomal = pclass.map(lambda s: "Ribosomal proteins" in s).fillna(False)
    # Immunoglobulin and TCR loci are the adaptive repertoire, not targets: in
    # myeloma and lymphoma lines the clonal Ig transcript is among the highest
    # in the cell, which put IGLC1 and IGLJ1 at the top of an earlier shortlist.
    is_immunoglobulin = pclass.map(
        lambda s: bool({"Immunoglobulin genes", "T-cell receptor genes"} & s)).fillna(False)

    subcell = D.load_subcellular().set_index("gene")
    main_loc = subcell["Main location"].reindex(annot.index).fillna("")
    add_loc = subcell["Additional location"].reindex(annot.index).fillna("")
    extracellular = subcell["Extracellular location"].reindex(annot.index).fillna("")
    at_pm = (main_loc + ";" + add_loc).str.contains("Plasma membrane", case=False)
    secreted_anywhere = annot["Secretome location"].fillna("").str.contains(
        "secreted", case=False)
    has_localisation = (main_loc.str.strip() != "") | (add_loc.str.strip() != "")

    organelle_only = main_loc.str.contains(ORGANELLE_TERMS, case=False, regex=True) & ~at_pm
    nuclear = main_loc.str.contains("Nucle", case=False)

    # Additive evidence score rather than a rule chain. Every single rule we
    # tried either admitted junk or rejected real drugs: HPA's antibody-based
    # localisation frequently places a surface receptor in the Golgi
    # (biosynthetic pool) without detecting plasma membrane, so treating
    # "organelle-only" as a hard disqualifier rejects BCMA, DLL3 and MSLN -
    # while admitting anything with a membrane topology lets in mitochondrial
    # and ER proteins. Weighing the evidence keeps both failure modes bounded
    # and makes the trade-off auditable (see `audit_surface_filter`).
    score = pd.Series(0, index=annot.index, dtype=int)
    score += 3 * is_cd_marker.astype(int)          # CD antigens are surface by definition
    score += 3 * is_gpcr.astype(int)               # seven-TM receptors
    score += 2 * (is_secreted | secreted_anywhere).astype(int)
    score += 2 * (extracellular != "").astype(int)  # measured extracellular signal
    score += 2 * at_pm.astype(int)                 # plasma-membrane localisation
    score += 1 * is_membrane.astype(int)           # membrane topology prediction
    # Penalties: no extracellular epitope, whatever the topology says.
    score -= 3 * is_tf.astype(int)
    score -= 3 * is_ribosomal.astype(int)
    score -= 6 * is_immunoglobulin.astype(int)
    score -= 2 * (nuclear & ~at_pm).astype(int)
    score -= 2 * (organelle_only & ~at_pm).astype(int)
    score -= 1 * is_intracellular.astype(int)
    # No localisation evidence at all is itself disqualifying: a gene whose
    # measured normal expression is near zero because it is poorly detected
    # looks identical to one that is genuinely tumour-restricted, and rewarding
    # that put PIGY, CYB561D2 and NRROS at the top of an earlier ranking.
    # ...but only when nothing else vouches for the protein. Secreted ligands
    # (VEGFA, ANGPT2) frequently have no immunofluorescence localisation at all,
    # and penalising that would reject the anti-angiogenic arm of real
    # bispecifics.
    strong_class = is_cd_marker | is_gpcr | is_secreted | secreted_anywhere
    score -= 2 * (~has_localisation & ~strong_class).astype(int)

    confidence = pd.Series("none", index=annot.index, dtype=object)
    confidence[score >= 1] = "medium"
    confidence[score >= 3] = "high"

    return pd.DataFrame({
        "is_membrane": is_membrane,
        "is_secreted": is_secreted,
        "at_plasma_membrane": at_pm.fillna(False),
        "organelle_only": organelle_only.fillna(False),
        "has_localisation": has_localisation.fillna(False),
        "surface_evidence_score": score,
        "surface_confidence": confidence,
        "surface_high": (confidence == "high"),
    })


def _expression_summary(genes: pd.Index) -> pd.DataFrame:
    mean_t, frac = D.tumour_by_cancer_type()
    normal = D.load_normal_sc().reindex(genes).fillna(0.0)
    mean_t = mean_t.reindex(genes).fillna(0.0)
    frac = frac.reindex(genes).fillna(0.0)

    vital = [c for c in C.VITAL_CELL_TYPES if c in normal.columns]
    tum_max = np.log1p(mean_t.max(axis=1))
    norm_max = np.log1p(normal.max(axis=1))
    norm_vital = np.log1p(normal[vital].max(axis=1))

    return pd.DataFrame({
        "best_tumour_line_fraction": frac.max(axis=1),
        "best_cancer_type": frac.idxmax(axis=1),
        "n_cancer_types_hit": (frac >= C.MIN_TUMOUR_LINE_FRACTION).sum(axis=1),
        "tumour_max_log": tum_max,
        "normal_max_log": norm_max,
        "normal_vital_max_log": norm_vital,
        "normal_breadth": (normal >= C.EXPR_EXPRESSED).sum(axis=1) / normal.shape[1],
        "window_log": tum_max - norm_max,
        "window_vital_log": tum_max - norm_vital,
    })


def build_universe(verbose: bool = True) -> pd.DataFrame:
    annot = D.load_annotations().set_index("gene")
    acc = _accessibility(annot)
    # Reachable means graded high or medium; "none" (including anything with
    # organelle-only localisation) is not a credible antibody target.
    accessible = acc.surface_confidence.isin(["high", "medium"])
    summ = _expression_summary(annot.index)

    eligible = accessible & (summ.best_tumour_line_fraction >= C.TUMOUR_ELIGIBILITY_FRACTION)

    anchors = (EFFECTOR_GENES | CHECKPOINT_GENES) & set(annot.index)
    known = known_pairs_df()
    known_genes = (set(known.gene_a) | set(known.gene_b)) & set(annot.index)

    genes = sorted(set(annot.index[eligible]) | anchors | known_genes)

    uni = pd.DataFrame(index=pd.Index(genes, name="gene")).join(acc).join(summ)
    uni["role"] = "TAA"
    uni.loc[uni.index.isin(EFFECTOR_GENES), "role"] = "effector"
    uni.loc[uni.index.isin(CHECKPOINT_GENES), "role"] = "checkpoint"
    uni["in_known_bsab"] = uni.index.isin(known_genes)
    uni["eligible_on_expression"] = uni.index.isin(annot.index[eligible])
    uni = uni.reset_index()

    if verbose:
        n_taa = int((uni.role == "TAA").sum())
        print(f"universe: {len(uni)} genes ({n_taa} TAA, "
              f"{int((uni.role == 'effector').sum())} effector, "
              f"{int((uni.role == 'checkpoint').sum())} checkpoint)")
        print("  surface confidence: "
              + ", ".join(f"{k} {v}" for k, v in uni.surface_confidence.value_counts().items()))
        forced = uni[uni.in_known_bsab & ~uni.eligible_on_expression]
        if len(forced):
            print(f"  retained despite low cell-line expression: {list(forced.gene)}")
    return uni


def audit_surface_filter(verbose: bool = True) -> dict:
    """How many validated bispecific targets the surface filter would reject.

    Reported rather than tuned away: it is the filter's false-negative rate on
    the only ground truth available.
    """
    annot = D.load_annotations().set_index("gene")
    acc = _accessibility(annot)
    known = known_pairs_df()
    genes = sorted((set(known.gene_a) | set(known.gene_b)) & set(annot.index))
    conf = acc.surface_confidence.reindex(genes)
    out = {
        "known_genes": len(genes),
        "high": int((conf == "high").sum()),
        "medium": int((conf == "medium").sum()),
        "rejected": int((conf == "none").sum()),
        "rejected_genes": sorted(conf.index[conf == "none"]),
        "high_only_would_miss": sorted(conf.index[conf != "high"]),
    }
    if verbose:
        print(f"surface filter on {out['known_genes']} validated targets: "
              f"{out['high']} high, {out['medium']} medium, {out['rejected']} rejected")
        if out["rejected_genes"]:
            print("  rejected:", out["rejected_genes"])
    return out


def benchmark_priors(universe: pd.DataFrame) -> pd.DataFrame:
    """Do hand-built target-quality priors beat chance at finding real targets?

    Each prior ranks the tumour antigens in the universe; we ask what fraction
    of the antigens used in documented bispecifics land in the top-k. The
    chance baseline is k / n, so `enrichment` is the ratio to chance. Values
    near 1 mean the prior carries essentially no information.
    """
    taa = universe[universe.role == "TAA"].set_index("gene")
    n = len(taa)
    known = taa.index[taa.in_known_bsab]
    if len(known) == 0:
        return pd.DataFrame()

    priors = {
        "prevalence": 2.0 * taa.best_tumour_line_fraction + 0.3 * np.log1p(taa.n_cancer_types_hit),
        "low_normal_expression": -taa.normal_max_log - 1.5 * taa.normal_breadth,
        "therapeutic_window": taa.window_log + 0.5 * taa.window_vital_log,
        "window_plus_prevalence": (1.5 * taa.best_tumour_line_fraction + 0.5 * taa.window_log
                                   - 0.35 * taa.normal_vital_max_log - 0.3 * taa.normal_breadth),
    }

    rows = []
    for name, s in priors.items():
        order = s.sort_values(ascending=False).index
        rank = {g: i for i, g in enumerate(order)}
        rec: dict[str, object] = {"prior": name}
        for k in (300, 600, 1200):
            hit = sum(1 for g in known if rank.get(g, n) < k)
            recall = hit / len(known)
            rec[f"recall@{k}"] = round(recall, 3)
            rec[f"enrich@{k}"] = round(recall / (k / n), 2)
        rec["median_rank"] = int(np.median([rank.get(g, n) for g in known]))
        rec["n_taa"] = n
        rec["n_known"] = int(len(known))
        rows.append(rec)
    return pd.DataFrame(rows)


def build_pair_space(universe: pd.DataFrame, verbose: bool = True) -> pd.DataFrame:
    """Structurally bounded pair space (see module docstring)."""
    info = universe.set_index("gene")
    taa = info.index[info.role == "TAA"]
    anchors = info.index[info.role.isin(["effector", "checkpoint"])]

    top_taa = (info.loc[taa, "best_tumour_line_fraction"]
               .sort_values(ascending=False)
               .index[:C.N_TOP_TAA_FOR_TAA_PAIRS])
    # Antigens used in documented bispecifics stay eligible for TAA x TAA pairs
    # so that positives are not silently removed from that block.
    top_taa = sorted(set(top_taa) | set(info.index[info.in_known_bsab & (info.role == "TAA")]))

    blocks: list[pd.DataFrame] = []

    def add(a_list, b_list, block: str, same: bool):
        a = np.asarray(sorted(a_list), dtype=object)
        b = np.asarray(sorted(b_list), dtype=object)
        if same:
            ia, ib = np.triu_indices(len(a), k=1)
            left, right = a[ia], a[ib]
        else:
            left = np.repeat(a, len(b))
            right = np.tile(b, len(a))
        # Order each pair alphabetically so `pair` is a canonical key. numpy's
        # minimum/maximum have no loop for mixed-width unicode, hence np.where.
        swap = left > right
        lo = np.where(swap, right, left)
        hi = np.where(swap, left, right)
        keep = lo != hi
        blocks.append(pd.DataFrame({"gene_a": lo[keep], "gene_b": hi[keep], "block": block}))

    add(taa, anchors, "TAA x anchor", same=False)
    add(anchors, anchors, "anchor x anchor", same=True)
    add(top_taa, top_taa, "TAA x TAA", same=True)

    pairs = pd.concat(blocks, ignore_index=True)
    pairs["pair"] = pairs.gene_a.str.cat(pairs.gene_b, sep="|")
    # Drop effector x effector and de-duplicate across blocks.
    role = info.role
    ra = role.reindex(pairs.gene_a).to_numpy()
    rb = role.reindex(pairs.gene_b).to_numpy()
    pairs = pairs[~((ra == "effector") & (rb == "effector"))]
    pairs = pairs.drop_duplicates("pair").reset_index(drop=True)

    if verbose:
        print(f"pair space: {len(pairs):,} pairs")
        print(pairs.block.value_counts().to_string())
    return pairs
