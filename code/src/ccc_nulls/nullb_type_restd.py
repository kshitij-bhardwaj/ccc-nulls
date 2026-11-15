"""Null B candidate T2: restandardised cell-type-ID permutation.

Same H0 and the same permutation as T1 (cell-type identities as whole units), but the
p-value does not rely on ranking alone, which caps T1 at 1/K^2. Following Efron &
Tibshirani's restandardisation, the centre and spread of the null come from the
cell-type permutation: the mean and SD of the interaction's score over its *other*
tested (source, target) pairs. The shape comes from the central limit theorem, because a
cell-type mean of many cells is close to normal (the approximation FastCCC uses for
clusters of more than 30 cells). The leave-one-out statistic

    t = (s(x, y) - mean_others) / (sd_others * sqrt(1 + 1/m)),   m = number of other pairs

is referred to t_{m-1}, the prediction-interval test for one more draw from the same
distribution. The observed pair is left out of its own reference, so a large effect is not
capped by the sample size. Abundance cancels, because centre and spread are per
interaction.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .type_perm import KEY, loo_t_pvalues, tested_groups


def p_values(res: pd.DataFrame, score: str = "score_obs", min_ref: int = 3) -> np.ndarray:
    """T2 p-value per row of ``res`` (NaN if untested or < ``min_ref`` other tested pairs)."""
    p = np.full(len(res), np.nan)
    for rows, g in tested_groups(res, KEY):
        p[rows] = loo_t_pvalues(g[score].to_numpy(float), min_ref)
    return p
