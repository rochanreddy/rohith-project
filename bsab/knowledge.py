"""Curated knowledge base of bispecific-antibody (bsAb) target pairs.

These pairs are the *positive* labels for supervised / pairwise learning.
They are drawn from approved oncology bsAbs and well-documented clinical-stage
programs (public sources: FDA/EMA labels, ClinicalTrials.gov, company pipelines).
Gene symbols follow HGNC (e.g. CD3 -> CD3E, BCMA -> TNFRSF17, 4-1BB -> TNFRSF9).

Every pair not in this list is *unlabelled*, not a confirmed negative, so the
learning problem is positive-unlabelled (PU). The modelling code accounts for that.
"""

from __future__ import annotations

import pandas as pd

# (target_a, target_b, status, mechanism_class, example_molecule)
KNOWN_PAIRS: list[tuple[str, str, str, str, str]] = [
    # ---------------- Approved oncology bsAbs ----------------
    ("CD19", "CD3E", "approved", "T-cell engager", "blinatumomab"),
    ("EPCAM", "CD3E", "approved", "T-cell engager", "catumaxomab"),
    ("TNFRSF17", "CD3E", "approved", "T-cell engager", "teclistamab / elranatamab / linvoseltamab"),
    ("MS4A1", "CD3E", "approved", "T-cell engager", "mosunetuzumab / glofitamab / epcoritamab"),
    ("GPRC5D", "CD3E", "approved", "T-cell engager", "talquetamab"),
    ("DLL3", "CD3E", "approved", "T-cell engager", "tarlatamab"),
    ("PMEL", "CD3E", "approved", "T-cell engager (TCR)", "tebentafusp"),
    ("EGFR", "MET", "approved", "dual signalling blockade", "amivantamab"),
    ("ERBB2", "ERBB3", "approved", "dual signalling blockade", "zenocutuzumab"),
    ("PDCD1", "CTLA4", "approved", "dual checkpoint", "cadonilimab"),
    ("PDCD1", "VEGFA", "approved", "checkpoint + anti-angiogenic", "ivonescimab"),
    # ---------------- Clinical-stage T-cell / NK-cell engagers ----------------
    ("FCRL5", "CD3E", "clinical", "T-cell engager", "cevostamab"),
    ("CD33", "CD3E", "clinical", "T-cell engager", "AMG 330"),
    ("IL3RA", "CD3E", "clinical", "T-cell engager", "flotetuzumab / vibecotamab"),
    ("CLEC12A", "CD3E", "clinical", "T-cell engager", "tepoditamab"),
    ("FLT3", "CD3E", "clinical", "T-cell engager", "AMG 427"),
    ("CD38", "CD3E", "clinical", "T-cell engager", "AMG 424"),
    ("CD22", "CD3E", "clinical", "T-cell engager", "JNJ-75348780"),
    ("CEACAM5", "CD3E", "clinical", "T-cell engager", "cibisatamab"),
    ("MUC16", "CD3E", "clinical", "T-cell engager", "ubamatamab"),
    ("FOLH1", "CD3E", "clinical", "T-cell engager", "pasotuxizumab / acapatamab"),
    ("STEAP1", "CD3E", "clinical", "T-cell engager", "xaluritamig"),
    ("KLK2", "CD3E", "clinical", "T-cell engager", "pasritamig"),
    ("GPC3", "CD3E", "clinical", "T-cell engager", "ERY974"),
    ("SSTR2", "CD3E", "clinical", "T-cell engager", "tidutamab"),
    ("CDH3", "CD3E", "clinical", "T-cell engager", "PF-06671008"),
    ("GUCY2C", "CD3E", "clinical", "T-cell engager", "PF-07062119"),
    ("CLDN18", "CD3E", "clinical", "T-cell engager", "ASP2138"),
    ("CLDN6", "CD3E", "clinical", "T-cell engager", "CLDN6xCD3 bsAb"),
    ("ERBB2", "CD3E", "clinical", "T-cell engager", "runimotamab"),
    ("CD276", "CD3E", "clinical", "T-cell engager", "MGD009"),
    ("ROR1", "CD3E", "clinical", "T-cell engager", "NVG-111"),
    ("MSLN", "CD3E", "clinical", "T-cell engager", "MSLN-directed TCE"),
    ("PSCA", "CD3E", "clinical", "T-cell engager", "GEM3PSCA"),
    ("TNFRSF8", "FCGR3A", "clinical", "NK-cell engager", "AFM13"),
    # ---------------- Co-stimulatory (signal 2) bsAbs ----------------
    ("FOLH1", "CD28", "clinical", "co-stimulatory engager", "REGN5678"),
    ("MUC16", "CD28", "clinical", "co-stimulatory engager", "REGN5668"),
    ("EGFR", "CD28", "clinical", "co-stimulatory engager", "REGN7075"),
    ("CEACAM5", "TNFRSF9", "clinical", "tumour-targeted 4-1BB agonist", "CEA-4-1BBL"),
    ("FAP", "TNFRSF9", "clinical", "tumour-targeted 4-1BB agonist", "RG7827"),
    ("ERBB2", "TNFRSF9", "clinical", "tumour-targeted 4-1BB agonist", "cinrebafusp alfa"),
    ("TPBG", "TNFRSF9", "clinical", "tumour-targeted 4-1BB agonist", "ALG.APV-527"),
    ("CLDN18", "TNFRSF9", "clinical", "tumour-targeted 4-1BB agonist", "givastomig"),
    ("CD274", "TNFRSF9", "clinical", "checkpoint + 4-1BB agonist", "acasunlimab"),
    ("FAP", "CD40", "clinical", "tumour-targeted CD40 agonist", "RO7300490"),
    ("MSLN", "CD40", "clinical", "tumour-targeted CD40 agonist", "ABBV-428"),
    ("FAP", "TNFRSF10B", "clinical", "tumour-targeted DR5 agonist", "RG7386"),
    ("CTLA4", "TNFRSF4", "clinical", "checkpoint + OX40 agonist", "ATOR-1015"),
    # ---------------- Dual checkpoint / immunomodulatory ----------------
    ("CD274", "CTLA4", "clinical", "dual checkpoint", "erfonrilimab (KN046)"),
    ("PDCD1", "LAG3", "clinical", "dual checkpoint", "tebotelimab"),
    ("CD274", "LAG3", "clinical", "dual checkpoint", "FS118"),
    ("PDCD1", "HAVCR2", "clinical", "dual checkpoint", "lomvastomig"),
    ("PDCD1", "TIGIT", "clinical", "dual checkpoint", "rilvegostomig"),
    ("CD274", "CD47", "clinical", "checkpoint + phagocytosis", "IBI322"),
    ("CD274", "VEGFA", "clinical", "checkpoint + anti-angiogenic", "PM8002"),
    ("CD274", "TGFB1", "clinical", "checkpoint + TGF-beta trap", "bintrafusp alfa"),
    ("MS4A1", "CD47", "clinical", "tumour-targeted CD47 blockade", "IMM0306"),
    ("CD19", "CD47", "clinical", "tumour-targeted CD47 blockade", "TG-1801"),
    # ---------------- Dual tumour / signalling / angiogenesis ----------------
    ("EGFR", "ERBB3", "clinical", "dual signalling blockade", "izalontamab (SI-B001)"),
    ("EGFR", "LGR5", "clinical", "dual signalling blockade", "petosemtamab"),
    ("ERBB3", "IGF1R", "clinical", "dual signalling blockade", "istiratumab"),
    ("DLL4", "VEGFA", "clinical", "dual anti-angiogenic", "navicixizumab"),
    ("ANGPT2", "VEGFA", "clinical", "dual anti-angiogenic", "vanucizumab"),
    ("CD19", "CD22", "clinical", "dual tumour antigen (escape prevention)", "dual CD19/CD22 constructs"),
]


def known_pairs_df() -> pd.DataFrame:
    df = pd.DataFrame(KNOWN_PAIRS, columns=["gene_a", "gene_b", "status", "mechanism", "example"])
    df["pair"] = [pair_key(a, b) for a, b in zip(df.gene_a, df.gene_b)]
    return df.drop_duplicates("pair").reset_index(drop=True)


def pair_key(a: str, b: str) -> str:
    """Order-invariant identifier for a target pair."""
    return "|".join(sorted((a, b)))


# Genes whose primary role is to recruit/activate immune effector cells.
EFFECTOR_GENES = {"CD3E", "CD28", "TNFRSF9", "TNFRSF4", "CD40", "FCGR3A", "CD2", "CD27", "ICOS", "KLRK1", "NCR1"}
# Immune checkpoint / immunomodulatory genes.
CHECKPOINT_GENES = {"PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "CD47", "VSIR", "BTLA", "SIRPA", "PDCD1LG2", "NT5E", "CD96"}
