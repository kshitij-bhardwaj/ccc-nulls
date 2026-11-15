"""Multiple-testing and p-value combination."""
from __future__ import annotations

import numpy as np


def bh(p) -> np.ndarray:
    """Benjamini–Hochberg q-values."""
    p = np.asarray(p, float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    q = p[order] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(q, 1.0)
    return out


def cauchy_combine(pvals, weights=None) -> np.ndarray:
    """Cauchy combination test (Liu & Xie, JASA 2020).

    ``pvals`` is (n_tests, n_components). The combined p-value is valid in the tail
    under arbitrary dependence between components:
    ``T = sum_i w_i tan((0.5 - p_i) pi)``,  ``p = 1/2 - arctan(T) / pi``.
    """
    P = np.clip(np.atleast_2d(np.asarray(pvals, float)), 1e-15, 1 - 1e-15)
    w = np.full(P.shape[1], 1 / P.shape[1]) if weights is None else np.asarray(weights, float)
    if not np.isclose(w.sum(), 1):
        raise ValueError("weights must sum to 1")
    T = (w * np.tan((0.5 - P) * np.pi)).sum(1)
    return 0.5 - np.arctan(T) / np.pi
