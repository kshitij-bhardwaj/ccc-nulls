"""Null B (B1, empirical null): known-answer checks of the fit and the output contract."""
import numpy as np
import pytest

from ccc_nulls import schema
from tools.cellchat import adapter, nulls
from tools.cellchat.empirical_null import empirical_null, fit_null_mle


def test_mle_recovers_a_known_null():
    """Pure N(2, 1.5^2), n = 20,000. The MLE is unbiased; its precision depends on the width of the
    central interval: 20 samples give sd(sigma0) ~0.18 on the central 50% (pct0 = 0.25) and ~0.03 on
    the central 80% (pct0 = 0.10). The interval is a sensitivity parameter (README)."""
    fits = [fit_null_mle(np.random.default_rng(s).normal(2, 1.5, 20000), pct0=0.10) for s in range(5)]
    assert all(abs(f.delta0 - 2) < 0.1 and abs(f.sigma0 - 1.5) < 0.1 and f.pi0 > 0.9 for f in fits)
    wide = [fit_null_mle(np.random.default_rng(s).normal(2, 1.5, 20000), pct0=0.25).sigma0 for s in range(10)]
    assert abs(np.mean(wide) - 1.5) < 0.15


@pytest.mark.parametrize("method", ["mle", "cm"])
def test_empirical_null_controls_fdr_where_the_theoretical_null_fails(method):
    """90% null N(2, 1.5^2) + 10% signal at N(8, 1.5^2): q <= 0.05 calls keep the false-discovery
    proportion near 5% with the empirical null; the theoretical N(0, 1) null does not."""
    rng = np.random.default_rng(1)
    fdp_emp, fdp_theo = [], []
    for _ in range(5):
        z = np.r_[rng.normal(2, 1.5, 9000), rng.normal(8, 1.5, 1000)]
        is_null = np.r_[np.ones(9000, bool), np.zeros(1000, bool)]
        for res, store in [(empirical_null(z, method, with_lfdr=False), fdp_emp),
                           (empirical_null(z, "theoretical", with_lfdr=False), fdp_theo)]:
            call = res["q"] <= 0.05
            store.append((call & is_null).sum() / max(call.sum(), 1))
    assert np.mean(fdp_emp) < 0.08
    assert np.mean(fdp_theo) > 0.5


def test_null_b_output_contract(small_adata):
    inputs = adapter.prepare_inputs(small_adata, "label")
    df = adapter.native_to_schema(adapter.run_native(inputs, nboot=2, seed=1))
    out = nulls.null_b(small_adata, df, "label", B=50, seed=0)
    schema.validate_null(df, out, "b") if hasattr(schema, "validate_null") else None
    for col in ("p_null_b", "q_null_b", "lfdr_null_b"):
        v = out[col].to_numpy()
        assert np.all((v > 0) & (v <= 1)), col
    assert out.attrs["null_b"]["n_fit"] > 0


def test_randomized_p_is_uniform_with_ties():
    """Observed and B permuted values exchangeable, 60% of them tied at 0 (CellChat's point mass):
    the randomized p is U(0, 1); the usual (1 + #>=) / (B + 1) is not."""
    from scipy import stats as sstats

    from tools.cellchat.tail_pvalues import randomized_p

    rng = np.random.default_rng(0)
    B, n = 200, 4000
    v = rng.gamma(2, size=(n, B + 1)) * (rng.uniform(size=(n, B + 1)) > 0.6)
    obs, perm = v[:, 0], v[:, 1:]
    p = randomized_p((perm > obs[:, None]).sum(1), (perm == obs[:, None]).sum(1), B, rng.uniform(size=n))
    assert sstats.kstest(p, "uniform").pvalue > 0.01
    assert sstats.kstest((1 + (perm >= obs[:, None]).sum(1)) / (B + 1), "uniform").pvalue < 1e-6


@pytest.mark.parametrize("positive", ["gamma", "lognormal"])
def test_tail_models_recover_small_p(positive):
    """Permutation distribution = 70% zeros + a continuous positive part with a known tail. From
    B = 1,000 draws, estimate p = 1e-4 and 1e-5 (unreachable with 1,000 permutations); log10 error,
    median over 30 repeats. GPD (xi >= 0) stays within ~half an order of magnitude whatever the
    tail. Pearson III is exact for the gamma part but, misspecified for the lognormal one, it
    understates p = 1e-5 by more than 10x: the moment model is only as good as its family."""
    from scipy import stats as sstats

    from tools.cellchat.tail_pvalues import N_EXC, ad_critical_values, gpd_tail_p, pearson3_tail_p

    dist = sstats.gamma(2) if positive == "gamma" else sstats.lognorm(0.5)
    B, zero = 1000, 0.7
    table = ad_critical_values()
    for p_true in (1e-4, 1e-5):
        x0 = dist.isf(p_true / (1 - zero))
        err = {"gpd": [], "pearson3": []}
        for s in range(30):
            rng = np.random.default_rng(s)
            perm = dist.rvs(size=B, random_state=rng) * (rng.uniform(size=B) > zero)
            pos = perm[perm > 0]
            top = -np.sort(-perm)[:N_EXC]
            err["gpd"].append(np.log10(gpd_tail_p(x0, top, B, table) / p_true))
            err["pearson3"].append(np.log10(pearson3_tail_p(x0, len(pos), pos.mean(), pos.std(ddof=1),
                                                            sstats.skew(pos), B) / p_true))
        med = {m: np.median(e) for m, e in err.items()}
        assert -0.6 < med["gpd"] < 0.8, (p_true, med)
        if positive == "gamma":
            assert abs(med["pearson3"]) < 0.3, (p_true, med)
        elif p_true == 1e-5:
            assert med["pearson3"] < -1, (p_true, med)
