import numpy as np
import pandas as pd
import pytest
from scipy import sparse

import ccc_nulls as cn
from ccc_nulls._matrix import cluster_means
from ccc_nulls.simulate import candidates_from_pairs, label_truth, simulate


@pytest.fixture(scope="module")
def s1():
    sim = simulate(seed=0, n_genes=3000, n_pairs=150, structured=True, planted=False)
    cands = candidates_from_pairs(sim.X, sim.genes, sim.labels, sim.pairs)
    return sim, cands


# ------------------------------------------------------------------ stats
def test_bh_matches_hand_computation():
    p = np.array([0.01, 0.04, 0.03, 0.5])
    assert np.allclose(cn.bh(p), [0.04, 0.16 / 3, 0.16 / 3, 0.5])


def test_cauchy_identities():
    p = np.array([[0.02, 0.02], [0.3, 0.3]])
    assert np.allclose(cn.cauchy_combine(p), [0.02, 0.3])             # equal p -> same p
    assert np.allclose(cn.cauchy_combine(np.array([[0.07]])), [0.07])  # single component
    with pytest.raises(ValueError):
        cn.cauchy_combine(p, weights=[0.7, 0.7])


# ------------------------------------------------------------------ null moments
def test_perm_moments_match_monte_carlo():
    rng = np.random.default_rng(1)
    X = rng.gamma(2.0, 1.0, (120, 15))
    labels = np.repeat([0, 1, 2], [60, 40, 20])
    E, V = cn.perm_moments(X, labels, 3)
    sims = np.stack([cluster_means(X, rng.permutation(labels), 3) for _ in range(4000)])
    assert np.allclose(sims.mean(0), E, atol=0.03)
    assert np.allclose(sims.var(0), V, rtol=0.15)


def test_gene_strata_shape_and_range(s1):
    sim, _ = s1
    st = cn.gene_strata(sim.X)
    assert st.shape == (sim.X.shape[1],) and st.min() == 0 and st.max() == 29


# ------------------------------------------------------------------ calibration
def test_structured_null_null_a_inflated_null_b_calibrated(s1):
    sim, cands = s1
    res = cn.significance_layer(sim.X, sim.genes, sim.labels, cands, B=300, resource_genes=set(sim.pairs.values.ravel()))
    fpr_a, fpr_b = (res.p_rand < 0.05).mean(), (res.p_strat < 0.05).mean()
    assert fpr_a > 0.15, fpr_a                 # uniform null mis-calibrated under structure
    assert fpr_b < 0.08, fpr_b                 # stratified null near nominal
    assert (res.q_strat < 0.05).sum() == 0


def test_fully_permuted_both_calibrated():
    sim = simulate(seed=2, n_genes=3000, n_pairs=150, structured=False, planted=False, shuffle_labels=True)
    cands = candidates_from_pairs(sim.X, sim.genes, sim.labels, sim.pairs)
    res = cn.significance_layer(sim.X, sim.genes, sim.labels, cands, B=300)
    assert 0.02 < (res.p_rand < 0.05).mean() < 0.08
    assert 0.02 < (res.p_strat < 0.05).mean() < 0.08


def test_planted_signal_recovered_by_null_b():
    sim = simulate(seed=3, n_genes=3000, n_pairs=150, n_planted=25, plant_lfc=2.0)
    cands = candidates_from_pairs(sim.X, sim.genes, sim.labels, sim.pairs)
    res = cn.significance_layer(sim.X, sim.genes, sim.labels, cands, B=300, resource_genes=set(sim.pairs.values.ravel()))
    truth = label_truth(res, sim)
    called = res.q_strat < 0.05
    assert called[truth == "true"].mean() > 0.6
    assert (called & (truth == "null")).sum() <= max(1, 0.1 * called.sum())


# ------------------------------------------------------------------ API behaviour
def test_sparse_matches_dense(s1):
    sim, cands = s1
    sub = cands.iloc[:300]
    a = cn.significance_layer(sim.X, sim.genes, sim.labels, sub, B=50, seed=7)
    b = cn.significance_layer(sparse.csr_matrix(sim.X), sim.genes, sim.labels, sub, B=50, seed=7)
    assert np.allclose(a.p_strat, b.p_strat) and np.allclose(a.p_rand, b.p_rand)


def test_unknown_genes_dropped_with_warning(s1):
    sim, cands = s1
    bad = pd.concat([cands.iloc[:20], pd.DataFrame([{"ligand": "NOPE", "receptor": sim.genes[0], "source": "T0", "target": "T1"}])])
    with pytest.warns(UserWarning, match="dropping 1"):
        res = cn.significance_layer(sim.X, sim.genes, sim.labels, bad, B=20)
    assert len(res) == 20


def test_missing_columns_raise(s1):
    sim, cands = s1
    with pytest.raises(ValueError, match="missing columns"):
        cn.significance_layer(sim.X, sim.genes, sim.labels, cands.drop(columns="target"), B=10)


def test_flags_and_summary():
    f = cn.flag_overlap([0.01, 0.01, 0.2, 0.2], [0.01, 0.2, 0.01, 0.2])
    assert list(f) == ["both", "A_only", "B_only", "none"]
    s = cn.overlap_summary(pd.DataFrame({"flag": f}))
    assert s.to_dict() == {"both": 1, "A_only": 1, "B_only": 1, "none": 1}


def test_aggregate_complexes_uses_max_p():
    res = pd.DataFrame({
        "ligand_complex": ["L1", "L1"], "receptor_complex": ["R1_R2", "R1_R2"],
        "source": ["A", "A"], "target": ["B", "B"],
        "p_rand": [0.001, 0.2], "p_strat": [0.01, 0.03], "score_obs": [2.0, 1.0],
    })
    agg = cn.aggregate_complexes(res)
    assert len(agg) == 1 and agg.p_rand.iloc[0] == 0.2 and agg.p_strat.iloc[0] == 0.03


def test_liana_adapter_schema():
    liana_res = pd.DataFrame({
        "ligand": ["A", "A"], "ligand_complex": ["A", "A"], "receptor": ["B", "C"], "receptor_complex": ["B_C", "B_C"],
        "source": ["x", "x"], "target": ["y", "y"], "lr_means": [1.0, 2.0], "cellphone_pvals": [0.0, 0.1],
    })
    out = cn.adapters.from_liana(liana_res, method="cellphonedb")
    assert list(out.columns) == ["ligand", "receptor", "source", "target", "tool_score", "ligand_complex", "receptor_complex"]
    with pytest.raises(KeyError):
        cn.adapters.from_liana(liana_res, tool_score="nope")


def test_binned_labels_diagnostic_equals_null_a(s1):
    sim, cands = s1
    sub = cands.iloc[:400]
    types = np.array(sorted(set(sim.labels)))
    lab = np.searchsorted(types, sim.labels)
    gi = pd.Index(sim.genes)
    li, ri = gi.get_indexer(sub.ligand), gi.get_indexer(sub.receptor)
    si, ti = np.searchsorted(types, sub.source), np.searchsorted(types, sub.target)
    _, pa = cn.null_a(sim.X, lab, len(types), li, ri, si, ti, B=300, rng=0)
    _, pl = cn.null_b_binned_labels(sim.X, lab, len(types), li, ri, si, ti, B=300, rng=1)
    assert np.corrcoef(np.log(pa), np.log(pl))[0, 1] > 0.95
