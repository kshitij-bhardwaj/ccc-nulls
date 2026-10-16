import numpy as np
import pytest
from scipy import sparse

from ccc_nulls import schema
from tools.cellchat import adapter, nulls
from tools.cellchat import cellchat_py as cc

pytest.importorskip("rdata")     # reads CellChatDB from the package's .rda the first time


def test_liana_reference_schema(small_adata):
    pytest.importorskip("liana")
    df = adapter.liana_reference(small_adata, "label", n_perms=10)
    assert len(df) > 0
    assert (df.score_name == adapter.SCORE_NAME).all()
    schema.validate(df)
    assert df.p_native.notna().all()


@pytest.mark.parametrize("seed, call, expected", [
    (1, lambda r: r.runif(3).round(7).tolist(), [0.2655087, 0.3721239, 0.5728534]),
    (123, lambda r: r.runif(1).round(7).tolist(), [0.2875775]),
    (1, lambda r: r.sample_int(10).tolist(), [9, 4, 7, 1, 2, 5, 3, 10, 6, 8]),
    (42, lambda r: r.sample_int(10).tolist(), [1, 5, 10, 8, 2, 4, 6, 9, 7, 3]),
    (123, lambda r: r.sample_int(10).tolist(), [3, 10, 2, 8, 6, 9, 1, 7, 5, 4]),
])
def test_r_generator_reproduces_r(seed, call, expected):
    """set.seed(seed) followed by runif / sample in R >= 3.6 gives these values."""
    assert call(cc.RRandom(seed)) == expected


def test_trimean_matches_sorted_reference():
    rng = np.random.default_rng(0)
    X = sparse.random(500, 30, density=0.3, random_state=1, data_rvs=lambda k: rng.gamma(1, 1, k)).tocsc()
    codes = rng.integers(0, 5, 500)
    fast = nulls.TriMean(X)(codes, 5)
    literal = cc.tri_mean_by_group(X.toarray().T, codes, 5).T
    assert np.array_equal(fast, literal)


@pytest.fixture(scope="module")
def inputs(small_adata):
    return adapter.prepare_inputs(small_adata, "label")


def test_native_pipeline_schema(inputs):
    df = adapter.native_to_schema(adapter.run_native(inputs, nboot=20, seed=1))
    schema.validate(df, require_p_native=True)
    assert (df.score > 0).any()


@pytest.mark.parametrize("version", ["2.2", "1.0"])
def test_null_a_reproduces_cellchat_exactly(small_adata, inputs, version):
    """With CellChat's seed, the vectorised Null A returns the literal port's prob and p bit for bit,
    for the current version (2.2) and the 2021 paper's (1.0)."""
    df = adapter.native_to_schema(adapter.run_native(inputs, nboot=20, seed=1, version=version))
    res = nulls.null_a(small_adata, df, "label", B=20, seed=1, rng="cellchat", version=version)
    assert np.array_equal(res.score_null_a.to_numpy(), df.score.to_numpy())
    assert np.array_equal(res.p_null_a_cellchat.to_numpy(), df.p_native.to_numpy())


def test_stats_and_collapse_quantiles():
    """stats::quantile interpolates (1 - h) a + h b only where b != a; collapse a + h (b - a)."""
    rng = np.random.default_rng(3)
    X = sparse.random(400, 20, density=0.4, random_state=2, data_rvs=lambda k: rng.gamma(1, 1, k)).tocsc()
    codes = rng.integers(0, 4, 400)
    for q in ("stats", "collapse"):
        assert np.array_equal(nulls.TriMean(X, quantile=q)(codes, 4), cc.tri_mean_by_group(X.toarray().T, codes, 4, q).T)
