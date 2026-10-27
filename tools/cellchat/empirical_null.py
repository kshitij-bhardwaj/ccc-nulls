"""Efron's empirical null for a vector of z-scores: the core of Null B (B1).

Why: when thousands of hypotheses are tested and most z-scores sit far from N(0, 1), the
theoretical null is valid but answers a question that is false almost everywhere. Efron's remedy
is to estimate the null distribution from the central bulk of the z's themselves and call only
what stands out from that bulk (Efron 2004; 2010, ch. 6). In Null B the z's are each interaction's
Null-A statistic, so "stands out" means "stronger than typical interactions", not "cell-type specific".

Steps and sources (defaults follow the R package locfdr, Efron, Turnbull & Narasimhan):

1. Null fit, f0 = N(delta0, sigma0^2) with null proportion pi0, from the central z's only
   (the "zero assumption": nearly every z in the centre is null; Efron 2004, sec. 4).
     "mle": maximum likelihood of a truncated normal on the central interval, pi0 from the
            count inside it (Efron 2007, "Size, power and false discovery rates", sec. 4;
            locfdr nulltype = 1, its default).
     "cm":  central matching, a quadratic fit to log f(z) near its mode (Efron 2004, sec. 4;
            locfdr nulltype = 2). Used as a sensitivity check.
   Central interval = [q(pct0), q(1 - pct0)], pct0 = 1/4 (locfdr's default).
2. Tail-area Fdr for the upper tail (Efron 2010, eq. 2.13; Storey 2002):
       Fdr(z) = pi0 * P0(Z >= z) / P(Z >= z),
   computed as pi0 x Benjamini-Hochberg q-values of p0 = 1 - Phi((z - delta0) / sigma0).
   This is the decision rule (call if q <= alpha): it needs only the null fit and the empirical
   tail proportions, with no density estimate.
3. Local fdr, fdr(z) = pi0 f0(z) / f(z) (Efron 2005; 2010, ch. 5), reported as a diagnostic.
   The mixture density f comes from Lindsey's method: Poisson regression of histogram counts on a
   natural cubic spline with 7 degrees of freedom plus an intercept (Efron & Tibshirani 1996;
   locfdr: 120 bins, df = 7). The pct most extreme z's on each side are left out of that fit
   (locfdr's `pct` argument); beyond them fdr is held at its value at the edge.

Only the upper tail is communication, so z below delta0 always gets fdr = 1 and is never called.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy import optimize
from scipy import stats as sstats

from ccc_nulls.stats import bh


@dataclass
class NullFit:
    delta0: float
    sigma0: float
    pi0: float
    method: str
    n: int
    interval: tuple[float, float]

    def to_dict(self):
        return asdict(self)


def central_interval(z, pct0=0.25):
    return float(np.quantile(z, pct0)), float(np.quantile(z, 1 - pct0))


def log_normal_mass(lo, hi):
    """log(Phi(hi) - Phi(lo)) for lo < hi, without cancellation: subtract in whichever tail is
    larger (log1p(-exp(.)) of the ratio). The naive difference underflows to 0 far from the mean,
    which would hand the likelihood a spurious bonus."""
    from scipy.special import log_ndtr

    if lo > 0:                                    # both in the upper tail: use survival functions
        return log_ndtr(-lo) + np.log1p(-np.exp(log_ndtr(-hi) - log_ndtr(-lo)))
    return log_ndtr(hi) + np.log1p(-np.exp(log_ndtr(lo) - log_ndtr(hi)))


def fit_null_mle(z: np.ndarray, pct0: float = 0.25) -> NullFit:
    """Truncated-normal maximum likelihood on the central interval [a, b] (Efron 2007, sec. 4).

    Every z in [a, b] is assumed null, so they are a sample from N(delta0, sigma0^2) truncated to
    [a, b]. Their log-likelihood is sum log phi((z - d) / s) - n_in log s - n_in log Q, with
    Q = Phi((b - d) / s) - Phi((a - d) / s), computed stably. Then pi0 = (n_in / n) / Q, capped at 1.
    The search is bounded (d within the interval widened by its length, s between 1/20 and 5
    interval lengths) so small samples cannot drift to s -> infinity, where the truncated
    likelihood is flat.
    """
    a, b = central_interval(z, pct0)
    zin = z[(z >= a) & (z <= b)]
    width = b - a
    if width <= 0:
        raise ValueError(f"central interval has zero width: over {1 - 2 * pct0:.0%} of the z's are equal")

    def nll(par):
        d, s = par[0], np.exp(par[1])
        return -(sstats.norm.logpdf(zin, d, s).sum() - len(zin) * log_normal_mass((a - d) / s, (b - d) / s))

    start = [float(np.median(zin)), float(np.log(width / 1.349))]       # IQR / 1.349 = sd for a normal
    bounds = [(a - width, b + width), (np.log(width / 20), np.log(5 * width))]
    res = optimize.minimize(nll, start, method="L-BFGS-B", bounds=bounds)
    d, s = float(res.x[0]), float(np.exp(res.x[1]))
    q = np.exp(log_normal_mass((a - d) / s, (b - d) / s))
    return NullFit(d, s, min(1.0, (len(zin) / len(z)) / q), "mle", int(len(z)), (a, b))


def mixture_density(z: np.ndarray, nbins: int = 120, df: int = 7, pct: float = 0.005):
    """Lindsey's method (Efron & Tibshirani 1996): histogram the z's, regress the counts on a
    natural cubic spline by Poisson GLM, and divide by n x bin width. Returns f and its range."""
    import patsy
    import statsmodels.api as sm

    lo, hi = np.quantile(z, [pct, 1 - pct]) if pct > 0 else (np.min(z), np.max(z))
    inside = z[(z >= lo) & (z <= hi)]
    edges = np.linspace(lo, hi, nbins + 1)
    counts, _ = np.histogram(inside, edges)
    mids = (edges[:-1] + edges[1:]) / 2
    # cr(x, df + 1) without a separate intercept spans the same space as locfdr's ns(x, df) + intercept
    basis = patsy.dmatrix(f"cr(x, df={df + 1}) - 1", {"x": mids}, return_type="dataframe")
    glm = sm.GLM(counts, basis, family=sm.families.Poisson()).fit()
    info, beta, scale = basis.design_info, glm.params.to_numpy(), len(z) * (edges[1] - edges[0])

    def f(x):
        X = patsy.build_design_matrices([info], {"x": np.clip(np.atleast_1d(x), lo, hi)})[0]
        return np.exp(np.asarray(X) @ beta) / scale

    return f, (float(lo), float(hi))


def fit_null_cm(z: np.ndarray, pct0: float = 0.25) -> NullFit:
    """Central matching (Efron 2004, sec. 4): log f(z) ~ c0 + c1 z + c2 z^2 over the central
    interval matches log(pi0 phi((z - d) / s) / s), so s = (-2 c2)^(-1/2), d = -c1 / (2 c2) and
    pi0 = exp(c0 - c1^2 / (4 c2)) sqrt(2 pi) s."""
    f, _ = mixture_density(z)
    a, b = central_interval(z, pct0)
    grid = np.linspace(a, b, 50)
    c2, c1, c0 = np.polyfit(grid, np.log(f(grid)), 2)
    if c2 >= 0:
        raise ValueError("central matching failed: log f(z) is not concave on the central interval")
    s = float(np.sqrt(-1 / (2 * c2)))
    pi0 = float(min(1.0, np.exp(c0 - c1 ** 2 / (4 * c2)) * np.sqrt(2 * np.pi) * s))
    return NullFit(float(-c1 / (2 * c2)), s, pi0, "cm", int(len(z)), (a, b))


def theoretical_null(z: np.ndarray) -> NullFit:
    """N(0, 1) with pi0 = 1: what the permutation z would follow if Null A's hypothesis held.
    Used as the ablation that removes the empirical step."""
    return NullFit(0.0, 1.0, 1.0, "theoretical", int(len(z)), (float("nan"), float("nan")))


def empirical_null(z: np.ndarray, method: str = "mle", pct0: float = 0.25, with_lfdr: bool = True) -> dict:
    """Fit the null to the finite z's; return per-z p0, q (tail-area Fdr), lfdr and the fit.
    Non-finite z (no permutation variation, so nothing to test) get p = q = lfdr = 1."""
    z = np.asarray(z, float)
    ok = np.isfinite(z)
    zf = z[ok]
    fit = {"mle": fit_null_mle, "cm": fit_null_cm}[method](zf, pct0) if method != "theoretical" else theoretical_null(zf)
    p0 = np.maximum(sstats.norm.sf(zf, fit.delta0, fit.sigma0), 1e-300)
    out = dict(fit=fit, p=np.ones(len(z)), q=np.ones(len(z)), lfdr=np.ones(len(z)))
    out["p"][ok] = p0
    out["q"][ok] = np.minimum(1.0, fit.pi0 * bh(p0))
    if with_lfdr:
        f, (lo, hi) = mixture_density(zf)
        zc = np.clip(zf, lo, hi)                                          # hold fdr at the edge beyond the fit
        fdr = np.minimum(1.0, fit.pi0 * sstats.norm.pdf(zc, fit.delta0, fit.sigma0) / np.maximum(f(zc), 1e-300))
        out["lfdr"][ok] = np.where(zf < fit.delta0, 1.0, fdr)
    return out
