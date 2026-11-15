import shutil

import numpy as np
import pandas as pd
import pytest

from benchmarks.liana_nulls import power, spatial
from benchmarks.liana_nulls.stages import kept_vs_dropped
from benchmarks.liana_nulls.tools import CellPhoneDB, complex_key


def test_cellphonedb_permutations_are_cellphonedbs_shuffles():
    """np.random.seed(s), then one shuffle of the cell-type column per iteration (cellphonedb 5.0.1
    shuffle_meta, threads = 1): indices from our replay give the same labels."""
    labels = np.repeat(np.array(list("abcde")), 7)
    P = CellPhoneDB("unused.zip", B=4).permutations(len(labels), 4)
    np.random.seed(0)
    for b in range(4):
        expected = labels.copy()
        np.random.shuffle(expected)
        assert np.array_equal(labels[P[:, b]], expected)


def test_complex_key_sorts_subunits():
    assert complex_key(["ITGB1_ITGA2", "TGFB1"]).tolist() == ["ITGA2_ITGB1", "TGFB1"]


def test_label_truth_true_one_sided_ambiguous_null():
    plan = pd.DataFrame(dict(ligand=["L1"], receptor=["R1"], source=["a"], target=["b"]))
    rows = pd.DataFrame(dict(source=["a", "a", "c", "a", "a"], target=["b", "c", "b", "b", "b"],
                             ligand_complex=["L1", "L1", "L1", "L1", "X"],
                             receptor_complex=["R1", "R1", "R1", "R1_R9", "Y"]))
    assert power.label_truth(rows, plan).tolist() == ["true", "one_sided", "one_sided", "ambiguous", "null"]


def test_metrics_power_and_fdp():
    rows = pd.DataFrame(dict(call=[True, True, False, True]))
    truth = np.array(["true", "one_sided", "true", "ambiguous"])
    m = power.metrics(rows, truth, np.array([True, False, True, False]), 2, {"x": "call"})
    assert m["x"]["power"] == 0.5 and m["x"]["fdp"] == 0.5 and m["x"]["one_sided_false_calls"] == 1


def test_raise_gene_only_touches_its_cells():
    from scipy import sparse
    rng = np.random.default_rng(0)
    X = sparse.csc_matrix(rng.gamma(1, 1, (50, 3)) * (rng.random((50, 3)) < 0.5))
    before = X.toarray()
    power.raise_gene(X, 1, np.arange(10), 3.0, rng)
    after = X.toarray()
    assert np.array_equal(after[10:], before[10:]) and np.array_equal(after[:, [0, 2]], before[:, [0, 2]])
    assert after[:10, 1].sum() > before[:10, 1].sum()


def test_colocalisation_flags_correlated_types():
    rng = np.random.default_rng(1)
    a = rng.random(500)
    props = np.c_[a, a + 0.05 * rng.random(500), rng.random(500), rng.random(500), rng.random(500)]
    t = spatial.colocalisation(props / props.sum(1, keepdims=True), np.array(list("abcde")))
    z = t.set_index(["source", "target"]).z
    assert z[("a", "b")] == z.max() and np.isnan(z[("a", "a")])


def test_kept_vs_dropped_odds_ratio():
    rows = pd.DataFrame(dict(source=["a"] * 4 + ["c"] * 4, target=["b"] * 4 + ["d"] * 4,
                             call_null_a=True, call_null_b=[True, True, True, False, True, False, False, False]))
    truth = pd.DataFrame(dict(source=["a", "c"], target=["b", "d"], colocalised=[True, False]))
    r = kept_vs_dropped(rows, truth, n_boot=50)
    assert r["kept_colocalised"] == 3 and r["dropped_colocalised"] == 1 and r["odds_ratio"] == pytest.approx(9.0)


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="needs R with CellChat 2.1.2")
def test_r_cellchat_equals_the_port(small_adata):
    pytest.importorskip("rdata")
    proc = __import__("subprocess").run(["Rscript", "-e", "library(CellChat)"], capture_output=True)
    if proc.returncode:
        pytest.skip("CellChat is not installed in R")
    from tools.cellchat import adapter, native_r
    inputs = adapter.prepare_inputs(small_adata, "label")
    py = adapter.run_native(inputs, nboot=20, seed=1)
    py["features_sig"] = inputs["features_sig"]
    fid = native_r.fidelity(native_r.run_r(small_adata, "label", nboot=20, seed=1), py)
    assert fid["lr_identical"] and fid["prob_identical"] == fid["n_rows"] and fid["pval_identical"] == fid["n_rows"]


def test_reference_profiles_pick_each_types_markers():
    import anndata as ad
    from scipy import sparse
    X = np.zeros((40, 6), np.float32)
    X[:20, :3], X[20:, 3:] = 2.0, 2.0
    a = ad.AnnData(sparse.csr_matrix(X), obs=pd.DataFrame(dict(t=["a"] * 20 + ["b"] * 20), index=[str(i) for i in range(40)]))
    a.var_names = [f"g{i}" for i in range(6)]
    prof = spatial.reference_profiles(a, "t", n_markers=3)
    assert prof.shape == (6, 2) and prof.loc["g0", "a"] == pytest.approx(np.expm1(2.0)) and prof.loc["g0", "b"] == 0


def test_presto_tie_rank_is_single_precision():
    """presto averages a tie group in float32: 22,207 tied zeros rank 11103.630859375, not 11104."""
    from tools.cellchat.cellchat_py import _presto_tie_rank
    assert _presto_tie_rank(0, 22207) == 11103.630859375
    assert _presto_tie_rank(0, 101) == 51.0 and _presto_tie_rank(41, 1) == 42.0


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="needs R with presto")
def test_wilcoxauc_equals_presto_on_sparse_tied_genes(tmp_path):
    """Sparse genes (many tied zeros, tied non-zeros): our p-values are presto's. The ranks and U are
    presto's exactly; the last bit can differ through R's pnorm vs scipy's ndtr (<= 1e-16 here)."""
    import subprocess
    from scipy import sparse
    from tools.cellchat.cellchat_py import wilcoxauc
    rng = np.random.default_rng(3)
    N, G, K = 20000, 6, 4
    X = (rng.integers(0, 4, (N, G)) * (rng.random((N, G)) < [[0.0004, 0.002, 0.01, 0.1, 0.5, 0.9]])).astype(np.float64)
    codes = rng.integers(0, K, N)
    np.savetxt(tmp_path / "x.txt", X), np.savetxt(tmp_path / "g.txt", codes, fmt="%d")
    r = subprocess.run(["Rscript", "-e", f"X <- t(as.matrix(read.table('{tmp_path}/x.txt'))); y <- factor(scan('{tmp_path}/g.txt', quiet=TRUE));"
                        "rownames(X) <- paste0('g', 1:nrow(X)); w <- presto::wilcoxauc(X, y);"
                        "cat(sprintf('%.17g', w$pval), sep='\\n')"], capture_output=True, text=True)
    if r.returncode:
        pytest.skip("presto is not installed in R")
    p_r = np.array(r.stdout.split(), float).reshape(K, G)          # presto: features within groups
    p_py = wilcoxauc(sparse.csc_matrix(X), codes, K)["pval"]
    assert np.abs(p_r - p_py).max() <= 1e-15 and np.array_equal(p_r < 0.05, p_py < 0.05)


def test_ulm_is_the_slope_t_value():
    from scipy import stats
    from benchmarks.liana_nulls.downstream import ulm
    rng = np.random.default_rng(0)
    w, Y = rng.normal(size=200), rng.normal(size=(200, 3))
    Y[:, 0] += 0.5 * w
    fit = stats.linregress(w, Y[:, 0])
    assert ulm(Y, w)[0] == pytest.approx(fit.slope / fit.stderr)


def test_lr_pathway_needs_every_subunit_with_one_sign():
    from benchmarks.liana_nulls.downstream import lr_pathways
    W = pd.DataFrame({"TGFb": [1.0, 2.0, 0.5, -1.0], "EGFR": [1.0, -1.0, 0.0, 1.0]}, index=["L1", "R1", "R2", "L2"])
    lr = pd.DataFrame(dict(ligand_complex=["L1", "L1", "L2"], receptor_complex=["R1", "R1_R2", "R1"]))
    m = lr_pathways(lr, W)
    assert set(map(tuple, m[["ligand_complex", "receptor_complex", "footprint", "sign"]].to_numpy().tolist())) == {
        ("L1", "R1", "TGFb", 1), ("L1", "R1_R2", "TGFb", 1)}          # mixed signs and zero weights are left out
