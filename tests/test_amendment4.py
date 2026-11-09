"""AMENDMENT_4: T2 switched on in every LIANA+ run of AMENDMENT_2's code."""
import warnings

import pytest

li = pytest.importorskip("liana")
pytest.importorskip("anndata")

from ccc_nulls import liana_consensus as lc  # noqa: E402

from test_liana_consensus import N_PERMS, SEED, adata  # noqa: E402,F401


@pytest.fixture
def switched_on():
    from benchmarks.liana_final import amendment2 as A2
    from benchmarks.liana_final import amendment4 as A4
    saved = (lc.run, lc.variants, lc.NULLS, A2.Results, A2.stability, A2.table, A2.ALL)
    A4.switch_on()
    yield A4
    lc.run, lc.variants, lc.NULLS, A2.Results, A2.stability, A2.table, A2.ALL = saved


def test_every_run_gets_t2(adata, switched_on):  # noqa: F811
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lrs, store = lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED, tools_only=True)
    assert all({"p_t2", "t2_tested"} <= set(d.columns) for d in store.values())


def test_null_b_consensus_also_builds_t2(adata, switched_on):  # noqa: F811
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lrs, store = lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED)
        rra = lc.variants(lrs, store, nulls=("null_a", "null_b"))
        mean = lc.variants(lrs, store, nulls=("null_b",), aggregate_method="mean")
    assert {"swap_t2", "extra_t2", "swap_null_b"} <= set(rra)
    assert "swap_t2" not in mean
