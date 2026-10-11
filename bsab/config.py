"""Paths and tunable constants for the bispecific target-pair pipeline."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = ROOT / "data" / "raw"
DATA_PROC = ROOT / "data" / "processed"
SCRNA_DIR = DATA_RAW / "scrna"
RESULTS = ROOT / "results"
FIGURES = RESULTS / "figures"

for _d in (DATA_PROC, RESULTS, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)

# ---- HPA source files ----
F_HPA_ANNOT = DATA_RAW / "proteinatlas.tsv"
F_HPA_SC = DATA_RAW / "rna_single_cell_type.tsv"
F_HPA_SUBCELL = DATA_RAW / "subcellular_location.tsv"
F_HPA_CELLLINE = DATA_RAW / "rna_celline.tsv"
F_CELLLINE_MAP = DATA_RAW / "celline_cancer_map.tsv"

# ---- Expression thresholds ----
# HPA treats nTPM/nCPM >= 1 as "detected"; >= 10 is a commonly used
# "clearly expressed" cut-off for cell-surface target triage.
EXPR_DETECTED = 1.0
EXPR_EXPRESSED = 10.0
EXPR_HIGH = 50.0

# Single-cell detection: a gene counts as "on" in a cell if raw count > 0.
SC_MIN_MALIGNANT_CELLS = 200  # skip datasets with fewer malignant cells

# ---- Candidate universe size control ----
# Universe membership: expressed (>= EXPR_EXPRESSED) in at least this fraction
# of the cell lines of at least one cancer type. Deliberately permissive - it
# is set low enough to retain every tumour antigen used in a documented
# bispecific (the lowest, ANGPT2, sits at 0.14).
TUMOUR_ELIGIBILITY_FRACTION = 0.10
# Stricter cut used for the "how many cancer types does this target hit" feature.
MIN_TUMOUR_LINE_FRACTION = 0.25
MIN_CANCER_TYPE_LINES = 6  # cancer types with fewer cell lines are dropped

# The pair space is bounded structurally rather than by a quality prior (which
# we measured to be near-chance at recovering validated targets). Every TAA is
# paired with every immune anchor; TAA x TAA pairs are limited to the most
# prevalent antigens by tumour-line positivity, a single transparent criterion.
N_TOP_TAA_FOR_TAA_PAIRS = 400

# ---- Normal cell types considered dose-limiting if targeted ----
# Damage to these lineages drives the classic on-target/off-tumour toxicities
# (cardiac, hepatic, neuro, pulmonary, renal, haematopoietic).
VITAL_CELL_TYPES = [
    "cardiomyocytes",
    "hepatocytes",
    "cholangiocytes",
    "brain excitatory neurons",
    "brain inhibitory neurons",
    "other brain neurons",
    "astrocytes",
    "oligodendrocytes",
    "alveolar cells type 1",
    "alveolar cells type 2",
    "proximal tubule cells",
    "distal convoluted tubule cells",
    "podocytes",
    "loop of henle epithelial cells",
    "hematopoietic stem cells",
    "erythrocyte progenitors",
    "megakaryocyte-erythroid progenitors",
    "monocyte progenitors",
    "neutrophil progenitors",
    "vascular endothelial cells",
    "pancreatic islet cells",
    "enterocytes",
    "colonocytes",
    "respiratory ciliated cells",
    "basal keratinocytes",
    "rod photoreceptor cells",
    "cone photoreceptor cells",
]

RANDOM_SEED = 20261011
