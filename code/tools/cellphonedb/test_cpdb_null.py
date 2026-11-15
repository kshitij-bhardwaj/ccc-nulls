"""Tests for the CellPhoneDB-faithful statistic and Null A (no database download needed)."""
import numpy as np
import pandas as pd
import pytest

from tools.cellphonedb import cpdb_null


def toy(seed=0, n_cells=60, signal=True):
    """Six proteins, one two-subunit complex, four interactions, three cell types."""
    rng = np.random.default_rng(seed)
    types = np.array(["B", "T", "Mono"])
    labels = np.repeat(types, n_cells // 3)
    rows = ["p1", "p2", "p3", "p4", "p5", "p6"]
    X = rng.gamma(0.6, 1.0, size=(len(rows), n_cells)) * (rng.random((len(rows), n_cells)) < 0.4)
    if signal:
        X[0, labels == "B"] += 2.0          # p1 high in B
        X[1, labels == "T"] += 2.0          # p2 high in T
    X[5] = 0.0                              # p6 never expressed: interactions with it score 0
    counts = pd.DataFrame(X.astype(np.float32), index=rows, columns=[f"c{i:03d}" for i in range(n_cells)])
    complexes = {"cx": [2, 3]}             # cx = p3 + p4
    interactions = pd.DataFrame({
        "multidata_1_id": ["p1", "cx", "p5", "p1"],
        "multidata_2_id": ["p2", "p1", "p6", "cx"],
        "id_cp_interaction": ["I1", "I2", "I3", "I4"],
    })
    cats = pd.Categorical(labels)
    data = cpdb_null.CpdbData(counts=counts.to_numpy(), row_ids=np.array(rows),
                              row_genes=np.array([r.upper() for r in rows]), codes=cats.codes,
                              cell_types=list(cats.categories), complexes=complexes, interactions=interactions)
    return data, counts, labels


def test_matches_cellphonedb_helpers():
    """Our means, fractions, scores and gate equal CellPhoneDB's own helper functions."""
    helper = pytest.importorskip("cellphonedb.src.core.methods.cpdb_statistical_analysis_helper")
    data, counts, labels = toy()
    meta = pd.DataFrame({"cell_type": labels}, index=counts.columns)
    clusters = helper.build_clusters(meta, counts, data.complexes, skip_percent=False)
    combos = helper.get_cluster_combinations(clusters["names"])
    inter = data.interactions[["multidata_1_id", "multidata_2_id"]]
    ref_mean = helper.mean_analysis(inter, clusters, combos, "|")
    ref_gate = helper.percent_analysis(clusters, 0.1, inter, combos, "|")

    means, fracs, _ = cpdb_null.feature_stats(data)
    feat = list(data.row_ids) + data.complex_ids
    assert list(clusters["names"]) == data.cell_types
    np.testing.assert_allclose(means, clusters["means"].loc[feat, data.cell_types].to_numpy(float), atol=1e-6)
    np.testing.assert_allclose(fracs, clusters["percents"].loc[feat, data.cell_types].to_numpy(float), atol=1e-12)
    assert list(ref_mean.columns) == data.pairs
    np.testing.assert_allclose(cpdb_null.pair_scores(means, data.i1, data.i2), ref_mean.to_numpy(), atol=1e-6)
    np.testing.assert_array_equal(cpdb_null.pair_gate(fracs, data.i1, data.i2, 0.1), ref_gate.to_numpy() == 1)


def test_statistic_rules_by_hand():
    data, counts, labels = toy()
    means, fracs, rep = cpdb_null.feature_stats(data)
    X = counts.to_numpy(float)
    for k, t in enumerate(data.cell_types):
        cells = labels == t
        sub = X[[2, 3]][:, cells].mean(1)
        assert means[len(data.row_ids), k] == pytest.approx(sub.min())                  # complex = min subunit
        assert data.row_genes[rep[len(data.row_ids), k]] == ["P3", "P4"][int(sub.argmin())]
        assert fracs[len(data.row_ids), k] == pytest.approx((X[[2, 3]][:, cells] > 0).mean(1).min())
    s = cpdb_null.pair_scores(means, data.i1, data.i2)
    assert (s[2] == 0).all()                                                            # p6 never expressed
    gate = cpdb_null.pair_gate(np.array([[0.1, 0.11]]), np.array([0]), np.array([0]), 0.1)
    assert gate.tolist() == [[False, False, False, True]]                               # strict '>'


def test_null_a_pvalues():
    data, _, _ = toy()
    res = cpdb_null.null_a(data, B=200, seed=1)
    assert res.p.shape == res.score.shape == (4, 9)
    assert np.isnan(res.p[~res.tested]).all() and not np.isnan(res.p[res.tested]).any()
    assert (res.p[res.tested] >= 1 / 201).all() and (res.p[res.tested] <= 1).all()
    assert not res.tested[2].any()                                                      # I3 scores 0 everywhere
    b_t = data.pairs.index("B|T")
    assert res.tested[0, b_t] and res.p[0, b_t] < 0.05                                  # planted p1(B) -> p2(T)
    cp = cpdb_null.null_a(data, B=200, seed=1, pvalue="cpdb")
    ok = res.tested
    np.testing.assert_array_less(cp.p[ok] - 1e-12, res.p[ok])                           # #>/B <= (1+#>=)/(B+1)
    again = cpdb_null.null_a(data, B=200, seed=1)
    np.testing.assert_array_equal(np.nan_to_num(again.p, nan=-1), np.nan_to_num(res.p, nan=-1))


def test_null_a_is_calibrated_without_signal():
    """With labels unrelated to expression, p-values of tested pairs are close to uniform."""
    ps = []
    for seed in range(6):
        data, _, _ = toy(seed=seed, n_cells=90, signal=False)
        res = cpdb_null.null_a(data, B=200, seed=100 + seed)
        ps.append(res.p[res.tested])
    p = np.concatenate(ps)
    assert len(p) > 50
    assert 0.0 <= (p <= 0.05).mean() <= 0.15
    assert 0.35 <= np.median(p) <= 0.65


def test_bad_inputs():
    data, counts, labels = toy()
    with pytest.raises(ValueError, match="separator"):
        cpdb_null.CpdbData(counts=counts.to_numpy(), row_ids=data.row_ids, row_genes=data.row_genes,
                           codes=data.codes, cell_types=["B|x", "T", "Mono"], complexes=data.complexes,
                           interactions=data.interactions)
    with pytest.raises(ValueError, match="pvalue"):
        cpdb_null.null_a(data, B=5, pvalue="exact")
