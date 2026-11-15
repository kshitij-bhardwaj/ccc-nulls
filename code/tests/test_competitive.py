"""Null B2 (ccc_nulls.competitive): the pieces against direct computations, and the correction on a
correlated null where the old decoy premise fails."""
import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats as sstats

from ccc_nulls import competitive as cp


def correlated_null(seed=1, N=3000, G=600, K=5, n_cand=100):
    """Random labels; decoys share a strong cell factor, candidates barely load on it."""
    rng = np.random.default_rng(seed)
    f = rng.normal(size=N)
    load = np.r_[np.full(n_cand, 0.05), np.full(G - n_cand, 0.6)]
    X = np.maximum(0, rng.uniform(0.3, 2, G) + rng.normal(size=(N, G)) * np.sqrt(1 - load**2) + f[:, None] * load)
    return X.astype(np.float32), rng.integers(0, K, N), [f"g{i}" for i in range(G)], K, n_cand


def test_welch_t_matches_scipy():
    X, codes, _, K, _ = correlated_null(N=800, G=40, n_cand=10)
    t, frac = cp.welch_t(sparse.csr_matrix(X), codes, K)
    for k in range(K):
        ref = sstats.ttest_ind(X[codes == k].astype(float), X[codes != k].astype(float), equal_var=False).statistic
        np.testing.assert_allclose(t[k], ref, rtol=1e-8)
        np.testing.assert_allclose(frac[k], (X[codes == k] > 0).mean(0))


def test_residuals_give_within_type_correlation():
    X, codes, _, K, _ = correlated_null(N=800, G=30, n_cand=10)
    U = cp.residuals(X, codes, K, n_cells=10_000)
    R = X.astype(np.float64).copy()
    for k in range(K):
        R[codes == k] -= R[codes == k].mean(0)
    np.testing.assert_allclose(U.T @ U, np.corrcoef(R.T), atol=1e-5)


def test_correction_restores_the_level_the_old_premise_loses():
    rates = []
    for rep in range(4):
        X, codes, genes, K, n = correlated_null(seed=10 + rep)
        gs = cp.gene_stats(X, genes, codes, K)
        groups = cp.decoy_groups(gs, genes[:n])
        idx = np.concatenate([np.arange(n)] + [np.asarray(g) for g in groups.values()])
        gs.with_residuals(X, codes, idx)
        sides = pd.DataFrame({"complex": genes[:n] * K, "type": np.repeat(np.arange(K), n)})
        st = cp.side_tests(gs, sides, groups)
        rates.append(((st.p < 0.05).mean(), (st.p_uncorrected < 0.05).mean()))
    corrected, uncorrected = np.mean(rates, 0)
    assert corrected <= 0.06, corrected
    assert uncorrected >= 0.08, uncorrected


def test_row_p_is_the_larger_side_p_and_gate_fails_to_one():
    X, codes, genes, K, n = correlated_null(seed=3, N=1500, G=400)
    gs = cp.gene_stats(X, genes, codes, K)
    rows = pd.DataFrame({"ligand_complex": ["g1", "g2_g3", "missing"], "receptor_complex": ["g4", "g5", "g6"],
                         "source": [0, 1, 2], "target": [1, 2, 3]})
    out = cp.null_b2(X, gs, codes, rows, genes[:n])
    np.testing.assert_array_equal(out.p_b2, np.maximum(out.p_L_b2, out.p_R_b2))
    assert out.p_b2.iat[2] == 1.0 and not out.gate_b2.iat[2]
