"""Tool scores, written against a types x genes matrix of cluster means ``M``.

Every score takes index arrays ``(li, ri, si, ti)`` (ligand gene, receptor gene,
source type, target type) and returns one value per candidate, so the same function
evaluates observed candidates, permuted data and decoy pairs.
"""
from __future__ import annotations

import numpy as np


def score_mean(M, li, ri, si, ti):
    """CellPhoneDB statistic: mean of ligand mean in source and receptor mean in target."""
    return 0.5 * (M[si, li] + M[ti, ri])


def score_product(M, li, ri, si, ti):
    """SingleCellSignalR / Connectome-style magnitude: product of the two means."""
    return M[si, li] * M[ti, ri]


SCORES = {"mean": score_mean, "product": score_product}


def get_score(score):
    if callable(score):
        return score
    try:
        return SCORES[score]
    except KeyError:
        raise ValueError(f"unknown score {score!r}; use one of {sorted(SCORES)} or a callable") from None


def standardized(f, M, E, V, li, ri, si, ti):
    """Tool score centred and scaled by its exact Null-A moments.

    ``E`` and ``V`` are the permutation mean / variance of every cluster mean
    (see :func:`ccc_nulls.nulls.perm_moments`). Ligand and receptor cluster means are
    treated as independent, which holds for distinct genes up to a small covariance
    through shared cells. Standardising puts genes of every abundance on one scale
    before they are ranked against decoys.
    """
    a, b = M[si, li], M[ti, ri]
    ea, eb, va, vb = E[si, li], E[ti, ri], V[si, li], V[ti, ri]
    if f is score_mean:
        return (0.5 * (a + b) - 0.5 * (ea + eb)) / np.sqrt(0.25 * (va + vb) + 1e-12)
    if f is score_product:
        return (a * b - ea * eb) / np.sqrt(ea**2 * vb + eb**2 * va + va * vb + 1e-12)
    raise ValueError("standardize=True supports the built-in 'mean' and 'product' scores; "
                     "pass standardize=False for a custom score")
