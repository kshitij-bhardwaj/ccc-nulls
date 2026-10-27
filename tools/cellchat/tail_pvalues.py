"""Permutation p-values for Null B (B1): exact where the B permutations resolve them, modelled beyond.

Why this step exists. Efron's empirical null (``empirical_null.py``) needs z-values that are N(0, 1)
for an interaction whose Null-A hypothesis is true. CellChat's prob is 0 in most permutations, so the
standardized score (prob - perm mean) / perm sd is not normal under the null: on fully shuffled
Tabula Sapiens it puts hundreds of interactions beyond z = 4 (README, "B1 calibration"). The
z-value that is exactly N(0, 1) under the null is the probit of an exact permutation p-value,
z = Phi^-1(1 - p), the construction Efron uses for every test he converts to z (Efron 2010, sec. 2.1).

Three p-values, all on CellChat's own prob and Null A's own permutations:

"permutation"  p = (n_greater + V (n_tie + 1)) / (B + 1), V ~ U(0, 1): the observed labelling is one
               of B + 1 exchangeable labellings, so its rank is uniform, and breaking ties at random
               makes p exactly U(0, 1) under the null (randomized test; Lehmann & Romano 2005,
               sec. 3.2 and 15.2). Ties matter here: an interaction with prob = 0 ties with every
               permutation that also gives 0. Its resolution ends at 1 / (B + 1).

Where fewer than M_MIN = 10 permutations reach the observed prob, the permutation p is too coarse
(Knijnenburg et al. 2009 show it is accurate from 10 exceedances on). There, either:

"gpd"          (option 1) Generalized Pareto tail (Knijnenburg, Wessels, Reinders & Shmulevich 2009,
               Bioinformatics 25:i161, "Fewer permutations, more accurate P-values"). By the
               Pickands-Balkema-de Haan theorem the excesses over a high threshold follow a GPD. Fit
               it by maximum likelihood to the N_exc = 250 largest permuted probs, threshold t halfway
               to the next one; check the fit with the Anderson-Darling test at 5%, and drop 10
               exceedances while it fails. Then p = (N_exc / B) * P_GPD(X > prob - t).
               One change, ours: the shape is held at xi >= 0. With a few hundred exceedances the MLE
               of xi is biased below 0 for exponential-type tails; xi < 0 means a finite upper
               endpoint, and extrapolating past it returns p = 0 (test_null_b.py: median 20x too small
               at p = 1e-4, often 0 at 1e-5). xi >= 0 assumes a tail no lighter than exponential, which
               overstates p when the true tail is lighter: a conservative error.
"pearson3"     (option 2) Moment-matched Pearson type III (a shifted gamma) on the first three
               permutation moments, the classical approximation to permutation distributions
               (Mielke & Berry 2007, "Permutation Methods", sec. 2.3; Zhou, Barry & Wright 2013,
               Biostatistics 14:97). CellChat's permutation distribution is a point mass at 0 plus a
               continuous part, so the moments are those of the positive part and
               p = (n_pos / B) * P_PIII(X >= prob | X > 0) (our two-part form; the zero mass is counted,
               not modelled).

Both models need at least MIN_TAIL = 30 positive permuted probs (a GPD fit on fewer is unstable, and
a third moment more so). Below that the conservative (n_greater + n_tie + 1) / (B + 1) is used
(Phipson & Smyth 2010), which is never below the randomized p, so the test stays valid.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import stats as sstats

M_MIN = 10          # Knijnenburg 2009: permutation p accurate from 10 exceedances on
N_EXC = 250         # their first number of exceedances, lowered in steps of STEP while the fit fails
STEP = 10
MIN_TAIL = 30       # ours: fewest positive permuted probs a tail model is fitted to
AD_LEVEL = 0.05     # Knijnenburg 2009: Anderson-Darling goodness of fit at 5%
AD_TABLE = Path(__file__).resolve().parent / "results" / "gpd_ad_critical.json"
P_FLOOR = 1e-300


def randomized_p(n_greater, n_tie, B, u):
    """Exact permutation p, U(0, 1) under the null: ties with the observed value broken at random."""
    return (n_greater + u * (n_tie + 1)) / (B + 1)


# ------------------------------------------------------------------ option 1: GPD tail
def anderson_darling(u):
    """A^2 of probability-transformed data u (Stephens 1974)."""
    u = np.clip(np.sort(u), 1e-12, 1 - 1e-12)
    n = len(u)
    i = np.arange(1, n + 1)
    return -n - np.mean((2 * i - 1) * (np.log(u) + np.log1p(-u[::-1])))


def ad_critical_values(n=100, sims=1000, seed=0):
    """5% critical value of A^2 for a GPD whose shape and scale are estimated, on a grid of shapes.

    The null distribution of A^2 with fitted parameters depends on the shape (Choulakian & Stephens
    2001, Technometrics 43:478, whose table this reproduces by parametric bootstrap). Computed once
    and cached in results/gpd_ad_critical.json.
    """
    if AD_TABLE.exists():
        t = json.loads(AD_TABLE.read_text())
        return np.array(t["shape"]), np.array(t["critical"])
    rng = np.random.default_rng(seed)
    shapes = np.round(np.arange(-0.5, 1.01, 0.1), 1)
    crit = []
    for xi in shapes:
        a2 = []
        for _ in range(sims):
            x = sstats.genpareto.rvs(xi, size=n, random_state=rng)
            c, _, sc = sstats.genpareto.fit(x, floc=0)
            a2.append(anderson_darling(sstats.genpareto.cdf(x, c, 0, sc)))
        crit.append(float(np.quantile(a2, 1 - AD_LEVEL)))
    AD_TABLE.write_text(json.dumps(dict(shape=shapes.tolist(), critical=crit, n=n, sims=sims, seed=seed,
                                        level=AD_LEVEL), indent=1))
    return shapes, np.array(crit)


def gpd_tail_p(x0, top, B, table):
    """Knijnenburg et al. (2009): p of the observed prob x0 from the descending largest permuted
    probs ``top`` (positive ones only). Returns NaN if no fit passes with >= MIN_TAIL exceedances."""
    shapes, crit = table
    y = np.append(top[top > 0], 0.0)                     # 0 below: threshold under the smallest positive
    n_exc = min(N_EXC, len(y) - 1)
    while n_exc >= MIN_TAIL:
        t = (y[n_exc - 1] + y[n_exc]) / 2                # halfway between the n_exc-th and next value
        exc = y[:n_exc] - t
        xi, _, sc = sstats.genpareto.fit(exc, floc=0)
        if xi < 0:                                       # exponential tail: MLE of the scale is the mean
            xi, sc = 0.0, float(exc.mean())
        a2 = anderson_darling(sstats.genpareto.cdf(exc, xi, 0, sc))
        if a2 <= np.interp(xi, shapes, crit):            # shapes outside the grid use its end values
            return max(P_FLOOR, n_exc / B * sstats.genpareto.sf(x0 - t, xi, 0, sc))
        n_exc -= STEP
    return np.nan


# ------------------------------------------------------------------ option 2: Pearson type III
def pearson3_tail_p(x0, n_pos, mean, sd, skew, B):
    """Two-part moment approximation: P(X >= x0) = P(X > 0) P_PIII(X >= x0 | X > 0), with the
    Pearson type III matched to the mean, sd and skewness of the positive permuted probs."""
    if n_pos < MIN_TAIL or sd <= 0:
        return np.nan
    return max(P_FLOOR, n_pos / B * sstats.pearson3.sf(x0, skew, loc=mean, scale=sd))


# ------------------------------------------------------------------ per interaction
def permutation_p(res, engine="gpd", seed=0):
    """p-value of every row of a ``nulls.null_a`` result. Returns (p, how), how in {"permutation",
    "tail", "capped"}: randomized permutation p, tail model, or conservative p where no model fits."""
    B = res.attrs["null_a"]["B"]
    prob = res["score_null_a"].to_numpy()
    k = res["n_exceed_null_a"].to_numpy()
    tie = res["n_tie_null_a"].to_numpy()
    p = randomized_p(k, tie, B, np.random.default_rng(seed).uniform(size=len(res)))
    how = np.full(len(res), "permutation", dtype=object)
    if engine == "permutation":
        return p, how
    tail = np.flatnonzero((prob > 0) & (k + tie < M_MIN))
    table = ad_critical_values() if engine == "gpd" else None
    for r in tail:
        if engine == "gpd":
            pt = gpd_tail_p(prob[r], res["perm_top_null_a"].iat[r], B, table)
        else:
            pt = pearson3_tail_p(prob[r], res["perm_npos_null_a"].iat[r], res["perm_pos_mean_null_a"].iat[r],
                                 res["perm_pos_sd_null_a"].iat[r], res["perm_pos_skew_null_a"].iat[r], B)
        if np.isfinite(pt):
            p[r], how[r] = pt, "tail"
        else:
            p[r], how[r] = (k[r] + tie[r] + 1) / (B + 1), "capped"
    return p, how
