"""Shared pieces for the cell-type-ID permutation nulls (Null B candidates).

These nulls permute cell-type identities as whole units instead of individual cells.
Each gene keeps its profile across cell types and every cell type keeps its size. Only
which type plays "source" and which plays "target" is randomised. Relabelling types
leaves a tool's score at (a, b) unchanged, so enumerating every relabelling is the same
as reading the tool's own score at every (source, target) pair. That makes the nulls
tool-agnostic: they need the tool's score for each interaction at every cell-type pair,
and nothing else.

Input everywhere: a long table (e.g. from tools/cellphonedb ``combine``) with one row per
interaction x (source, target), a score column and ``tested`` (the tool's expression gate).
Only tested rows enter the reference sets, so candidates and references face the same
filter.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

KEY = ["ligand_complex", "receptor_complex"]     # one interaction = one complex pair


def tested_groups(res: pd.DataFrame, by: list[str]):
    """Yield (row positions in ``res``, group) for each group of tested rows."""
    pos = pd.Series(np.arange(len(res)), index=res.index)
    for _, g in res.loc[res["tested"].astype(bool)].groupby(by, sort=False):
        yield pos[g.index].to_numpy(), g


def loo_t_pvalues(s: np.ndarray, min_ref: int = 3) -> np.ndarray:
    """Upper-tail p-value of each value against all the *other* values in ``s``.

    Leave-one-out studentised score ``t = (s_i - mean_-i) / (sd_-i * sqrt(1 + 1/m))`` with
    m = len(s) - 1, referred to a t distribution with m - 1 degrees of freedom (the
    prediction-interval test for a new draw from the same distribution). NaN when fewer than
    ``min_ref`` reference values. p = 1 when the references have no spread and s_i does not
    exceed them, and p = 0 when it does.
    """
    s = np.asarray(s, float)
    m = len(s) - 1
    if m < min_ref:
        return np.full(len(s), np.nan)
    mu = (s.sum() - s) / m
    var = ((s ** 2).sum() - s ** 2 - m * mu ** 2) / (m - 1)
    var = np.maximum(var, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (s - mu) / np.sqrt(var * (1 + 1 / m))
    p = stats.t.sf(t, m - 1)
    flat = var <= 1e-12 * np.maximum(mu ** 2, 1e-12)
    return np.where(flat, np.where(s > mu, 0.0, 1.0), p)
