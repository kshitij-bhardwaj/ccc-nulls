"""AMENDMENT_3: T2 (PR #11) added to LIANA+'s captured tools, and the swap_t2 / extra_t2 consensus variants."""
import warnings

import numpy as np
import pytest

li = pytest.importorskip("liana")
pytest.importorskip("anndata")

from ccc_nulls import liana_consensus as lc  # noqa: E402
from ccc_nulls.nullb_type_restd import p_values  # noqa: E402

from test_liana_consensus import N_PERMS, SEED, adata  # noqa: E402,F401


@pytest.fixture(scope="module")
def t2_run(adata):  # noqa: F811
    from benchmarks.liana_final import amendment3 as A3
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lrs, store = lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED)
    A3.add_t2(store)
    return lrs, store


def test_p_t2_is_t2_on_the_tool_magnitude(t2_run):
    _, store = t2_run
    for d in store.values():
        u = d.drop_duplicates(lc.KEY).reset_index(drop=True)
        ref = p_values(u.assign(tested=True), score="magnitude")
        got = u.merge(d.drop_duplicates(lc.KEY)[lc.KEY + ["p_t2", "t2_tested"]], on=lc.KEY, suffixes=("_u", ""))
        assert np.allclose(got.p_t2, np.where(np.isfinite(ref), ref, 1.0))
        assert (got.t2_tested == np.isfinite(ref)).all()
        assert d.p_t2.between(0, 1).all()


def test_t2_variants_leave_the_others_unchanged(t2_run):
    lrs, store = t2_run
    old = lc.variants(lrs, store, nulls=("null_a", "null_b"))
    new = lc.variants(lrs, store, nulls=("null_a", "null_b", "t2"))
    assert {"swap_t2", "extra_t2"} <= set(new)
    for v in old:
        a, b = old[v].sort_values(lc.KEY).reset_index(drop=True), new[v].sort_values(lc.KEY).reset_index(drop=True)
        assert a.equals(b), v
