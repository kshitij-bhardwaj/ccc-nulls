"""Null B (composite null, intersection-union test): the split score, the output contract and a
known-null check on real PBMC data."""
import numpy as np

from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat import cellchat_py as cc
from tools.cellchat.null_b_calibration import shuffle_genes, side_genes


def _candidates(adata):
    return adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(adata, "label"), nboot=2, seed=1))


def test_null_b_contract_and_same_score_as_null_a(small_adata):
    df = _candidates(small_adata)
    a = nulls.null_a(small_adata, df, "label", B=50, seed=3, rng="numpy")
    b = nulls.null_b(small_adata, df, "label", B=50, seed=3, rng="numpy", delta=(0.0, 0.05))
    assert np.array_equal(a.score_null_a.to_numpy(), b.score_null_b.to_numpy())     # the score is not changed
    pL, pR, p = b.p_L_null_b.to_numpy(), b.p_R_null_b.to_numpy(), b.p_null_b.to_numpy()
    assert np.array_equal(p, np.maximum(pL, pR))
    assert np.all(p[b.score_null_b.to_numpy() == 0] == 1)
    assert np.all(b["p_null_b_delta=0.05"].to_numpy() >= p)                          # interval null is stricter


def test_receptor_shuffle_is_a_known_null_for_null_b_not_for_null_a(small_adata):
    """Receptor-side genes moved to random cells: the ligand side keeps its cell-type signal, the
    receptor side has none, so "s signals to t" is false for every pair whose receptor side was
    shuffled. Null A still calls such pairs (the ligand alone drives the score); Null B calls none."""
    _, rec_only = side_genes(cc.subset_db(cc.load_db("human")))
    ad = shuffle_genes(small_adata, rec_only, np.random.default_rng(0))
    df = _candidates(ad)
    a = nulls.null_a(ad, df, "label", B=200, seed=1, rng="numpy")
    b = nulls.null_b(ad, df, "label", B=200, seed=1, rng="numpy")
    known_null = df.receptor_genes.map(lambda x: all(g in rec_only for g in nulls._genes(x))).to_numpy()
    family = nulls.testable(a, b)
    calls = {}
    for name, pv in [("A", a.p_null_a.to_numpy()), ("B", b.p_null_b.to_numpy())]:
        q = np.ones(len(df))
        q[family] = bh(pv[family])
        calls[name] = int(((q <= 0.05) & known_null).sum())
    assert calls["A"] > 0 and calls["B"] == 0, calls
