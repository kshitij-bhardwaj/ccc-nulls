"""ccc_nulls.composite: Null B as a composite null, on any tool's side-split score."""
import numpy as np
import pandas as pd
import pytest

from ccc_nulls import composite
from ccc_nulls._matrix import cluster_means
from ccc_nulls.nulls import null_a
from ccc_nulls.scores import score_mean, score_product


def _one_sided_data(seed=0, n=600, K=3):
    """Gene 0 (ligand) is high in type 0 only; gene 1 (receptor) has the same values in every type; gene 2
    (receptor) is high in type 1 only. So (0 -> 1 via genes 0, 1) has a specific sender and a
    non-specific receiver, a known null for "s signals to t"; (0 -> 1 via genes 0, 2) is a true edge."""
    rng = np.random.default_rng(seed)
    labels = np.repeat(np.arange(K), n // K)
    X = rng.gamma(2.0, 0.5, size=(n, 3))
    X[labels == 0, 0] += 3.0
    X[:, 1] = np.tile(X[: n // K, 1], K)
    X[labels == 1, 2] += 3.0
    return X, labels


@pytest.mark.parametrize("f", [score_mean, score_product])
def test_mean_side_score_equals_the_tool_score(f):
    X, labels = _one_sided_data()
    li, ri, si, ti = np.array([0, 0, 2]), np.array([1, 2, 1]), np.array([0, 0, 2]), np.array([1, 1, 0])
    side = composite.mean_side_score(X, li, ri, si, ti, f)
    M = side.summarize(labels, 3)
    assert np.allclose(side.score(M, M), f(cluster_means(X, labels, 3), li, ri, si, ti))


def test_one_specific_side_is_called_by_null_a_not_by_null_b():
    X, labels = _one_sided_data()
    li, ri, si, ti = np.array([0, 0]), np.array([1, 2]), np.array([0, 0]), np.array([1, 1])
    _, p_a = null_a(X, labels, 3, li, ri, si, ti, "mean", B=199, rng=1)
    res = composite.composite(composite.mean_side_score(X, li, ri, si, ti), labels, 3, B=199, seed=1)
    assert p_a[0] <= 0.01                         # the ligand alone drives Null A
    assert res["p"][0] > 0.5 and res["limiting_side"][0] == "R"
    assert res["p"][1] <= 0.01                     # both sides specific: Null B calls it
    assert np.array_equal(res["p"], np.maximum(res["p_L"], res["p_R"]))


def test_delta_ties_and_given_permutations():
    X, labels = _one_sided_data()
    side = composite.mean_side_score(X, np.array([0, 0]), np.array([1, 2]), np.array([0, 0]), np.array([1, 1]))
    perms = np.stack([np.random.default_rng(b).permutation(len(labels)) for b in range(50)], 1)
    r = composite.composite(side, labels, 3, permutations=perms, delta=(0.0, 0.5))
    assert r["B"] == 50 and r["deltas"] == [0.0, 0.5]
    assert np.all(r["p_delta"][0.5] >= r["p"])                                    # interval null is stricter
    r2 = composite.composite(side, labels, 3, permutations=perms, ties="ignore")
    assert np.all(r2["p"] <= r["p"])
    assert set(composite.to_columns(r)) >= {"score_null_b", "p_L_null_b", "p_R_null_b", "p_null_b",
                                            "limiting_side_null_b", "p_null_b_delta=0.5"}


def test_testable_is_any_labelling_with_a_positive_score():
    t = composite.testable([0, 0, 0, 1], [0, 2, 0, 0], [0, 0, 3, 0], [True, False, False, False])
    assert t.tolist() == [True, True, True, True]
    assert composite.testable([0, 0], [0, 0], [0, 0]).tolist() == [False, False]


def test_registry_lets_a_new_null_b_run_through_any_side_score():
    calls = []

    def every_row_null(side, codes, n_types, **kw):
        s = np.asarray(side.score(side.summarize(codes, n_types), side.summarize(codes, n_types)))
        calls.append(kw)
        return dict(score=s, p=np.ones(s.shape))

    assert composite.get("composite") is composite.composite
    assert composite.get(every_row_null) is every_row_null
    composite.register("test_every_row_null")(every_row_null)
    try:
        with pytest.raises(ValueError):
            composite.register("test_every_row_null")(every_row_null)
        X, labels = _one_sided_data()
        genes = ["L", "R1", "R2"]
        cand = pd.DataFrame(dict(ligand=["L", "L"], receptor=["R1", "R2"], source=["0", "0"], target=["1", "1"]))
        out = composite.null_b_table(X, genes, labels, cand, method="test_every_row_null", B=7)
        assert out.p_null_b.tolist() == [1.0, 1.0] and not out.call_null_b.any()
        assert calls[-1]["B"] == 7
    finally:
        composite.METHODS.pop("test_every_row_null")
    with pytest.raises(ValueError):
        composite.get("no_such_method")


def test_null_b_table_on_common_schema():
    X, labels = _one_sided_data()
    cand = pd.DataFrame(dict(ligand=["L", "L", "nope"], receptor=["R1", "R2", "R1"],
                             source=["0", "0", "0"], target=["1", "1", "1"], extra=[1, 2, 3]))
    with pytest.warns(UserWarning, match="dropping 1"):
        out = composite.null_b_table(X, ["L", "R1", "R2"], labels, cand, B=199, seed=2)
    assert len(out) == 2 and "extra" in out.columns
    assert out.call_null_b.tolist() == [False, True]
    assert out.attrs["null_b"]["method"] == "composite" and out.attrs["null_b"]["n_testable"] == 2
