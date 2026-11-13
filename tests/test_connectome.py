"""Connectome in the consensus (liana_consensus.capture_connectome, swap3_<n>): LIANA+'s outputs are unchanged,
and Connectome's specificity is replaced only in the swap3 variants."""
import warnings

import numpy as np
import pandas as pd
import pytest

li = pytest.importorskip("liana")
from ccc_nulls import liana_consensus as lc  # noqa: E402

from test_liana_consensus import N_PERMS, SEED, _merge, adata  # noqa: E402,F401  (the same toy data)


def _runs(a):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return (lc.run(a.copy(), "type", n_perms=N_PERMS, seed=SEED),
                lc.run(a.copy(), "type", n_perms=N_PERMS, seed=SEED, connectome=True))


@pytest.fixture(scope="module")
def both(adata):  # noqa: F811
    (lrs0, s0), (lrs1, s1) = _runs(adata)
    for s in (s0, s1):
        for d in s.values():
            d["p_b2"] = d["p_null_b"]          # B2 needs a decoy universe; not tested here
    return (lrs0, s0), (lrs1, s1)


def test_baseline_and_two_tool_variants_unchanged(both):
    (lrs0, s0), (lrs1, s1) = both
    V0, V1 = lc.variants(lrs0, s0, nulls=("null_a", "null_b")), lc.variants(lrs1, s1, nulls=("null_a", "null_b"))
    for v in V0:
        m = _merge(V0[v], V1[v])
        assert len(m) == len(V0[v]) == len(V1[v])
        np.testing.assert_array_equal(m.specificity_rank_x, m.specificity_rank_y)
        np.testing.assert_array_equal(m.magnitude_rank_x, m.magnitude_rank_y)
    for t in lc.TOOLS:
        pd.testing.assert_frame_equal(s0[t], s1[t])
    assert set(V1) - set(V0) == {"swap3_null_a", "swap3_null_b"}


def test_captured_connectome_is_liana_connectome(adata, both):  # noqa: F811
    _, (lrs1, s1) = both
    c = s1["Connectome"]
    ref = li.mt.connectome(adata.copy(), "type", use_raw=False, inplace=False)
    m = _merge(c, ref)
    assert len(m) == len(ref)
    np.testing.assert_allclose(m.magnitude, m.expr_prod)
    m2 = _merge(lrs1["Connectome"], ref)
    np.testing.assert_array_equal(m2.scaled_weight_x, m2.scaled_weight_y)   # LIANA+'s specificity, untouched
    assert ((c.p_null_a >= 0) & (c.p_null_a <= 1)).all()
    assert (c.p_null_b >= c[["p_L_null_b", "p_R_null_b"]].max(axis=1) - 1e-12).all()


def test_swap3_replaces_connectome_specificity(both):
    _, (lrs1, s1) = both
    V = lc.variants(lrs1, s1, nulls=("null_a", "null_b"))
    for n in ("null_a", "null_b"):
        cn = _merge(V[f"swap3_{n}"], s1["Connectome"])
        np.testing.assert_array_equal(cn.connectome_p, cn[f"p_{n}"])
    cp = _merge(V["swap3_null_a"], s1["CellPhoneDB"])
    np.testing.assert_array_equal(cp.cellphone_pvals, cp.p_null_a)        # null_a keeps the tools' own p
    cp = _merge(V["swap3_null_b"], s1["CellPhoneDB"])
    np.testing.assert_array_equal(cp.cellphone_pvals, cp.p_null_b)
