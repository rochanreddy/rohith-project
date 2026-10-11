"""Loaders for the Human Protein Atlas tables used by the pipeline.

Everything here is cached to parquet under data/processed/ because the raw
single-cell and cell-line tables are long-format and slow to re-pivot.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import config as C


def _cache(name: str, build, force: bool = False) -> pd.DataFrame:
    path = C.DATA_PROC / f"{name}.parquet"
    if path.exists() and not force:
        return pd.read_parquet(path)
    df = build()
    df.to_parquet(path)
    return df


# --------------------------------------------------------------------------
# Gene-level annotation
# --------------------------------------------------------------------------
def load_annotations(force: bool = False) -> pd.DataFrame:
    """One row per gene symbol with the annotation columns we use as features."""

    def build() -> pd.DataFrame:
        keep = [
            "Gene", "Ensembl", "Gene description", "Protein class", "Biological process",
            "Molecular function", "Disease involvement", "Subcellular location",
            "Subcellular main location", "Secretome location", "Secretome function",
            "RNA tissue specificity", "RNA tissue specificity score",
            "RNA single cell type specificity", "RNA single cell type specificity score",
            "RNA cancer specificity", "RNA cancer specificity score",
            "Interactions",
        ]
        df = pd.read_csv(C.F_HPA_ANNOT, sep="\t", low_memory=False, usecols=keep)
        df = df.rename(columns={"Gene": "gene"})
        # A handful of symbols appear twice (paralogue annotation); keep the
        # entry with the most populated fields.
        df["_fill"] = df.notna().sum(axis=1)
        df = df.sort_values("_fill", ascending=False).drop_duplicates("gene").drop(columns="_fill")
        return df.set_index("gene").sort_index().reset_index()

    return _cache("hpa_annotations", build, force)


def load_subcellular(force: bool = False) -> pd.DataFrame:
    def build() -> pd.DataFrame:
        df = pd.read_csv(C.F_HPA_SUBCELL, sep="\t", low_memory=False)
        df = df.rename(columns={"Gene name": "gene"})
        cols = ["gene", "Reliability", "Main location", "Additional location", "Extracellular location"]
        df = df[cols].drop_duplicates("gene")
        return df.reset_index(drop=True)

    return _cache("hpa_subcellular", build, force)


# --------------------------------------------------------------------------
# Normal single-cell-type expression (gene x normal cell type, nCPM)
# --------------------------------------------------------------------------
def load_normal_sc(force: bool = False) -> pd.DataFrame:
    """Wide matrix: index = gene symbol, columns = normal cell types, values nCPM."""

    def build() -> pd.DataFrame:
        df = pd.read_csv(
            C.F_HPA_SC, sep="\t",
            usecols=["Gene name", "Cell type", "nCPM"],
            dtype={"Gene name": "string", "Cell type": "category", "nCPM": "float32"},
        )
        wide = df.pivot_table(index="Gene name", columns="Cell type", values="nCPM",
                              aggfunc="max", observed=True)
        wide.index.name = "gene"
        wide.columns = [str(c) for c in wide.columns]
        return wide.astype("float32")

    return _cache("normal_sc_wide", build, force)


# --------------------------------------------------------------------------
# Cancer cell-line expression (gene x cell line, nTPM) + cancer-type map
# --------------------------------------------------------------------------
def load_celline_map() -> pd.DataFrame:
    df = pd.read_csv(C.F_CELLLINE_MAP, sep="\t", header=None,
                     names=["cell_line", "cancer_type"])
    df = df[~df.cancer_type.isin(["non-cancerous", "uncategorized"])]
    counts = df.cancer_type.value_counts()
    keep = counts[counts >= C.MIN_CANCER_TYPE_LINES].index
    return df[df.cancer_type.isin(keep)].drop_duplicates("cell_line").reset_index(drop=True)


def load_tumour_lines(force: bool = False) -> pd.DataFrame:
    """Wide matrix: index = gene symbol, columns = cancer cell line, values nTPM."""

    def build() -> pd.DataFrame:
        head = pd.read_csv(C.F_HPA_CELLLINE, sep="\t", nrows=5)
        gene_col = "Gene name" if "Gene name" in head.columns else head.columns[1]
        line_col = next(c for c in head.columns if "ell line" in c or c == "Cell line")
        # Prefer nTPM: HPA's between-sample normalised value, which is the one
        # comparable across cell lines. Plain TPM is not.
        val_col = next((c for c in ("nTPM", "pTPM", "TPM") if c in head.columns), None)
        if val_col is None:
            raise ValueError(f"no expression column in {C.F_HPA_CELLLINE}: {list(head.columns)}")
        df = pd.read_csv(
            C.F_HPA_CELLLINE, sep="\t",
            usecols=[gene_col, line_col, val_col],
            dtype={gene_col: "string", line_col: "category", val_col: "float32"},
        )
        wide = df.pivot_table(index=gene_col, columns=line_col, values=val_col,
                              aggfunc="max", observed=True)
        wide.index.name = "gene"
        wide.columns = [str(c) for c in wide.columns]
        return wide.astype("float32")

    return _cache("tumour_lines_wide", build, force)


def tumour_by_cancer_type(force: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-cancer-type summaries of tumour expression.

    Returns (mean_ntpm, frac_lines_expressed) with genes as rows and cancer
    types as columns. The fraction uses EXPR_EXPRESSED as the positivity cut-off,
    which matters more than the mean for antibody targeting: a target is useful
    when most tumour cells carry it, not when a few lines carry a lot.
    """
    lines = load_tumour_lines(force)
    cmap = load_celline_map()
    cmap = cmap[cmap.cell_line.isin(lines.columns)]
    means, fracs = {}, {}
    for ctype, grp in cmap.groupby("cancer_type"):
        sub = lines[list(grp.cell_line)]
        means[ctype] = sub.mean(axis=1)
        fracs[ctype] = (sub >= C.EXPR_EXPRESSED).mean(axis=1)
    return (pd.DataFrame(means).astype("float32"),
            pd.DataFrame(fracs).astype("float32"))


# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------
def tau_specificity(row: np.ndarray) -> float:
    """Yanai tau: 0 = ubiquitous, 1 = exclusive to one column."""
    x = np.asarray(row, dtype=float)
    x = x[~np.isnan(x)]
    if x.size < 2:
        return np.nan
    mx = x.max()
    if mx <= 0:
        return np.nan
    return float((1.0 - x / mx).sum() / (x.size - 1))


def split_terms(value) -> set[str]:
    if not isinstance(value, str) or not value.strip():
        return set()
    return {t.strip() for t in value.split(",") if t.strip()}
