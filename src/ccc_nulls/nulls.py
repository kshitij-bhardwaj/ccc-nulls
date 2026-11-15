"""The two null models.

Null A (random)      Permute cell-type labels uniformly and recompute the tool score.
                     This is the null CellPhoneDB and CellChat use.

Null B (stratified)  Keep the real labels, so cell-type sizes and every gene's profile
                     across cell types are preserved. Rank the candidate (L, R, X, Y)
                     against every decoy pair (L', R') whose genes share L's / R's
                     mean-expression decile x detection-rate stratum. A gene's overall
                     abundance then confers no advantage.

All functions take integer index arrays: ``li, ri`` (gene columns of ligand / receptor)
and ``si, ti`` (source / target cell-type codes). ``labels`` holds integer codes
0..n_types-1, one per cell.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ._matrix import cluster_means, columns, detection, detection_by_type, gene_mean_var
from .scores import get_score, standardized


# ------------------------------------------------------------------ strata
def gene_strata(X, n_mean_bins: int = 10, n_det_bins: int = 3) -> np.ndarray:
    """Stratum id per gene: mean-expression quantile bin x detection-rate bin.

    Detection bins are computed *within* each mean bin, so the strata are balanced.
    """
    mu, _ = gene_mean_var(X)
    det = detection(X)
    n_mean_bins = min(n_mean_bins, len(mu))
    mbin = pd.qcut(pd.Series(mu).rank(method="first"), n_mean_bins, labels=False).to_numpy()
    dbin = np.zeros_like(mbin)
    for b in range(n_mean_bins):
        idx = np.where(mbin == b)[0]
        if len(idx) >= n_det_bins:
            dbin[idx] = pd.qcut(pd.Series(det[idx]).rank(method="first"), n_det_bins, labels=False)
    return (mbin * n_det_bins + dbin).astype(np.int64)


# ------------------------------------------------------------------ Null A
def null_a(X, labels, n_types, li, ri, si, ti, score="mean", B=1000, rng=None):
    """Uniform label permutation.

    Returns ``(s_obs, p)`` with ``p = (1 + #{b : s(pi_b) >= s_obs}) / (B + 1)``.
    """
    rng = np.random.default_rng(rng)
    f = get_score(score)
    used, inv = np.unique(np.concatenate([li, ri]), return_inverse=True)
    li_, ri_ = inv[: len(li)], inv[len(li):]
    Xu = columns(X, used)                      # only genes that appear in candidates
    s_obs = f(cluster_means(Xu, labels, n_types), li_, ri_, si, ti)
    exceed = np.zeros(len(s_obs))
    for _ in range(B):
        M = cluster_means(Xu, rng.permutation(labels), n_types)
        exceed += f(M, li_, ri_, si, ti) >= s_obs
    return s_obs, (1 + exceed) / (B + 1)


def perm_moments(X, labels, n_types):
    """Exact mean and variance of every cluster mean under uniform label permutation.

    A permuted cluster of size n_k is a simple random sample without replacement, so
    ``Var = sigma^2 / n_k * (n - n_k) / (n - 1)`` (finite-population correction).
    Returns two types x genes arrays.
    """
    n = len(labels)
    nk = np.bincount(labels, minlength=n_types)[:, None].astype(float)
    mu, var = gene_mean_var(X)
    V = var[None, :] / nk * (n - nk) / (n - 1)
    return np.broadcast_to(mu[None, :], V.shape), V


# ------------------------------------------------------------------ Null B
def null_b(X, labels, n_types, li, ri, si, ti, score="mean", rng=None, strata=None,
           standardize=True, decoy_pool=None, min_frac=0.10, max_decoys=400,
           return_n_ref=False):
    """Expression-stratified null with exact enumeration over decoy pairs.

    Parameters
    ----------
    strata
        Stratum id per gene (default :func:`gene_strata`).
    standardize
        Rank the Null-A-standardised score instead of the raw score, so heterogeneity
        *within* a stratum does not cap power. Requires a built-in score.
    decoy_pool
        Boolean mask of genes allowed as decoys. Default: every gene that is not a
        candidate ligand or receptor. Pass ``~is_resource_gene`` to exclude the whole
        LR resource. LR genes may carry real signal and would contaminate the null.
    min_frac
        The expression filter the candidates passed. Decoys for source type X must be
        expressed in >= ``min_frac`` of X cells (likewise for the target), otherwise
        the filter itself biases candidates upward relative to decoys.
    max_decoys
        Cap on decoy genes per (stratum, cell type) group, subsampled at random.

    return_n_ref
        Also return the size of each candidate's reference set. ``1 / n_ref`` is the
        smallest p-value Null B can reach for that candidate.

    Returns
    -------
    (s_obs, p) or (s_obs, p, n_ref), where s_obs is the (standardised) score.

    Notes
    -----
    Under H0 the ligand is exchangeable with its group's decoys, and so is the receptor,
    so (L, R) is a uniform draw from (Lg + {L}) x (Rg + {R}). The p-value is the
    candidate's upper-tail rank in that whole product set, including the pairs that
    reuse L or R. Leaving those out lets one extreme ligand reach p ~ 1/N instead of
    ~ 1/|Lg|, which is anti-conservative.
    """
    rng = np.random.default_rng(rng)
    base = get_score(score)
    M = cluster_means(X, labels, n_types)
    if standardize:
        E, V = perm_moments(X, labels, n_types)

        def f(M_, l, r, s, t):
            return standardized(base, M_, E, V, l, r, s, t)
    else:
        f = base
    s_obs = f(M, li, ri, si, ti)

    strata = gene_strata(X) if strata is None else np.asarray(strata)
    if decoy_pool is None:
        decoy_pool = np.ones(X.shape[1], bool)
        decoy_pool[np.concatenate([li, ri])] = False
    frac = detection_by_type(X, labels, n_types)
    n_str = int(strata.max()) + 1
    G2 = n_str * n_types

    groups = []                                # index = stratum * n_types + type
    for st in range(n_str):
        in_st = strata == st
        for k in range(n_types):
            ok = in_st & (frac[k] >= min_frac)
            g = np.where(ok & decoy_pool)[0]
            if len(g) == 0:                    # fall back: stratum + filter, then stratum
                g = np.where(ok)[0]
            if len(g) == 0:
                g = np.where(in_st)[0]
            if max_decoys and len(g) > max_decoys:
                g = np.sort(rng.choice(g, max_decoys, replace=False))
            groups.append(g)

    lgrp = strata[li] * n_types + si
    rgrp = strata[ri] * n_types + ti
    keys = lgrp.astype(np.int64) * G2 + rgrp
    p = np.empty(len(li))
    n_ref = np.empty(len(li), dtype=np.int64)
    for key in np.unique(keys):
        idx = np.where(keys == key)[0]
        Lg, Rg = groups[key // G2], groups[key % G2]
        nL, nR, m = len(Lg), len(Rg), len(idx)
        s0, t0 = si[idx[0]], ti[idx[0]]
        obs = s_obs[idx]
        dd = np.sort(f(M, np.repeat(Lg, nR), np.tile(Rg, nL), np.full(nL * nR, s0), np.full(nL * nR, t0)))
        cnt = nL * nR - np.searchsorted(dd, obs, side="left")
        lr = f(M, np.repeat(li[idx], nR), np.tile(Rg, m), np.full(m * nR, s0), np.full(m * nR, t0)).reshape(m, nR)
        rl = f(M, np.tile(Lg, m), np.repeat(ri[idx], nL), np.full(m * nL, s0), np.full(m * nL, t0)).reshape(m, nL)
        cnt = cnt + (lr >= obs[:, None]).sum(1) + (rl >= obs[:, None]).sum(1)
        p[idx] = (1 + cnt) / ((nL + 1) * (nR + 1))
        n_ref[idx] = (nL + 1) * (nR + 1)
    return (s_obs, p, n_ref) if return_n_ref else (s_obs, p)


def null_b_binned_labels(X, labels, n_types, li, ri, si, ti, score="mean", B=1000, rng=None, strata=None):
    """Diagnostic: permute cell labels independently *within each gene stratum*.

    This is the literal reading of "permute labels within expression bins". Each gene
    still sees a uniform label permutation, so every candidate's marginal null equals
    Null A's, and the abundance bias is not removed. Kept to document that result
    (``tests/test_ccc_nulls.py::test_binned_labels_diagnostic_equals_null_a``).
    """
    rng = np.random.default_rng(rng)
    f = get_score(score)
    strata = gene_strata(X) if strata is None else np.asarray(strata)
    used, inv = np.unique(np.concatenate([li, ri]), return_inverse=True)
    li_, ri_ = inv[: len(li)], inv[len(li):]
    Xu, su = columns(X, used), strata[used]
    s_obs = f(cluster_means(Xu, labels, n_types), li_, ri_, si, ti)
    bins = [np.where(su == b)[0] for b in np.unique(su)]
    exceed = np.zeros(len(s_obs))
    M = np.empty((n_types, len(used)))
    for _ in range(B):
        for cols in bins:
            M[:, cols] = cluster_means(columns(Xu, cols), rng.permutation(labels), n_types)
        exceed += f(M, li_, ri_, si, ti) >= s_obs
    return s_obs, (1 + exceed) / (B + 1)
