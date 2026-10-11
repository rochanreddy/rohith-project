"""Chunked pairwise reductions over gene x condition expression matrices.

With ~5,600 genes a full gene-by-gene matrix is 127 MB per statistic and the
cell-line matrix (1,206 conditions) makes the naive loop hopeless. Every
statistic here is instead computed only for the pairs actually requested, in
row chunks, so memory stays flat in the number of pairs.
"""

from __future__ import annotations

import numpy as np

DEFAULT_CHUNK = 20_000


def _chunks(n: int, size: int):
    for start in range(0, n, size):
        yield start, min(start + size, n)


def max_of_min(P: np.ndarray, ia: np.ndarray, ib: np.ndarray,
               chunk: int = DEFAULT_CHUNK) -> np.ndarray:
    """For each pair, max over conditions of min(expr_a, expr_b).

    This is the AND-gate exposure: for a bispecific that needs both antigens on
    the same cell, what a condition contributes is bounded by the *lesser* of
    the two, and the risk is set by the worst such condition.
    """
    out = np.empty(len(ia), dtype=np.float32)
    for s, e in _chunks(len(ia), chunk):
        out[s:e] = np.minimum(P[ia[s:e]], P[ib[s:e]]).max(axis=1)
    return out


def min_of_max(P: np.ndarray, ia: np.ndarray, ib: np.ndarray,
               chunk: int = DEFAULT_CHUNK) -> np.ndarray:
    """For each pair, min over conditions of max(expr_a, expr_b).

    The OR-gate floor: coverage a bispecific retains if either antigen suffices.
    """
    out = np.empty(len(ia), dtype=np.float32)
    for s, e in _chunks(len(ia), chunk):
        out[s:e] = np.maximum(P[ia[s:e]], P[ib[s:e]]).min(axis=1)
    return out


def both_above(P: np.ndarray, ia: np.ndarray, ib: np.ndarray, thresh: float,
               chunk: int = DEFAULT_CHUNK) -> np.ndarray:
    """Count of conditions where both genes exceed `thresh`."""
    B = (P >= thresh).astype(np.float32)
    out = np.empty(len(ia), dtype=np.float32)
    for s, e in _chunks(len(ia), chunk):
        out[s:e] = (B[ia[s:e]] * B[ib[s:e]]).sum(axis=1)
    return out


def either_above(P: np.ndarray, ia: np.ndarray, ib: np.ndarray, thresh: float,
                 chunk: int = DEFAULT_CHUNK) -> np.ndarray:
    B = (P >= thresh).astype(np.float32)
    out = np.empty(len(ia), dtype=np.float32)
    for s, e in _chunks(len(ia), chunk):
        out[s:e] = np.maximum(B[ia[s:e]], B[ib[s:e]]).sum(axis=1)
    return out


def row_corr(P: np.ndarray, ia: np.ndarray, ib: np.ndarray,
             chunk: int = DEFAULT_CHUNK) -> np.ndarray:
    """Pearson correlation between the two genes' condition profiles.

    High correlation means the two antigens rise and fall together across
    conditions - informative but double-edged: it predicts co-expression on
    tumour, and also co-expression on whatever normal tissue carries them.
    """
    mu = P.mean(axis=1, keepdims=True)
    sd = P.std(axis=1, keepdims=True)
    Z = np.divide(P - mu, sd, out=np.zeros_like(P, dtype=np.float32),
                  where=sd > 1e-9)
    n = P.shape[1]
    out = np.empty(len(ia), dtype=np.float32)
    for s, e in _chunks(len(ia), chunk):
        out[s:e] = (Z[ia[s:e]] * Z[ib[s:e]]).sum(axis=1) / n
    return np.nan_to_num(out)
