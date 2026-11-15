"""NATMI in the consensus (liana_consensus.capture_natmi, swap4_<n>): LIANA+'s outputs are unchanged, every earlier
variant (Connectome's included) is the same with NATMI captured, NATMI's p-values are those of a direct computation
from the permutations, and NATMI's specificity is replaced only in the swap4 variants."""
import warnings

import numpy as np
import pandas as pd
import pytest

li = pytest.importorskip("liana")
from ccc_nulls import liana_consensus as lc  # noqa: E402

from test_liana_consensus import N_PERMS, SEED, _merge, adata  # noqa: E402,F401  (the same toy data)

NULLS = ("null_a", "null_b")


@pytest.fixture(scope="module")
def both(adata):  # noqa: F811
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r0 = lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED, connectome=True)
        r1 = lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED, connectome=True, natmi=True)
    for _, s in (r0, r1):
        for d in s.values():
            d["p_b2"] = d["p_null_b"]          # B2 needs a decoy universe; not tested here
    return r0, r1


def test_earlier_variants_and_tools_unchanged(both):
    (lrs0, s0), (lrs1, s1) = both
    V0, V1 = lc.variants(lrs0, s0, nulls=NULLS), lc.variants(lrs1, s1, nulls=NULLS)
    for v in V0:
        m = _merge(V0[v], V1[v])
        assert len(m) == len(V0[v]) == len(V1[v])
        np.testing.assert_array_equal(m.specificity_rank_x, m.specificity_rank_y)
        np.testing.assert_array_equal(m.magnitude_rank_x, m.magnitude_rank_y)
    for t in s0:
        pd.testing.assert_frame_equal(s0[t], s1[t])
    assert set(V1) - set(V0) == {"swap4_null_a", "swap4_null_b"}


def test_shared_permutations_change_nothing(adata):  # noqa: F811
    """Connectome with the shared cube equals Connectome with its own (both built from the same seed)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, s_shared = lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED, connectome=True)
        store = {}
        import liana as li_
        from liana.method.sc._rank_aggregate import AggregateClass, _rank_aggregate_meta
        agg = AggregateClass(_rank_aggregate_meta, methods=[lc.capture(li_.mt.cellphonedb, store), lc.capture_connectome(store)])
        with lc.lazy_permutations():
            agg(adata.copy(), groupby="type", resource_name="consensus", n_perms=N_PERMS, seed=SEED, expr_prop=0.1,
                consensus_opts=False, inplace=False, use_raw=False)
    for t in ("CellPhoneDB", "Connectome"):
        pd.testing.assert_frame_equal(s_shared[t], store[t])


def test_captured_natmi_is_liana_natmi(adata, both):  # noqa: F811
    _, (lrs1, s1) = both
    n = s1["NATMI"]
    ref = li.mt.natmi(adata.copy(), "type", use_raw=False, inplace=False)
    m = _merge(n, ref)
    assert len(m) == len(ref)
    np.testing.assert_allclose(m.expr_prod_x, m.expr_prod_y)
    ok = np.isfinite(m.spec_weight)
    np.testing.assert_allclose(m.magnitude[ok], m.spec_weight[ok], rtol=1e-6)   # spec_weight (LIANA+: float32)
    m2 = _merge(lrs1["NATMI"], ref)
    np.testing.assert_array_equal(m2.spec_weight_x.fillna(-1), m2.spec_weight_y.fillna(-1))   # LIANA+'s, untouched
    assert ((n.p_null_a >= 0) & (n.p_null_a <= 1)).all()
    np.testing.assert_array_equal(n.p_null_b, n[["p_L_null_b", "p_R_null_b"]].max(axis=1))


def test_natmi_p_from_the_permutations(adata, both):  # noqa: F811
    """NATMI's p-values recomputed directly: the cube LIANA+ builds (same seed), each row's types and genes."""
    _, (_, s1) = both
    import liana.method.sc._liana_pipe as pipe
    cubes = []
    inner = pipe._get_means_perms
    pipe._get_means_perms = lambda *a, **k: cubes.append(inner(*a, **k)) or cubes[-1]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED, natmi=True, tools_only=True)
    finally:
        pipe._get_means_perms = inner
    cube = cubes[0]
    types = list(pd.Categorical(adata.obs["type"]).categories)
    lab = adata.obs["type"].astype(str).to_numpy()
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    means = np.stack([X[lab == t].mean(0) for t in types])           # observed type means, as LIANA+'s
    tot_obs, tot = means.sum(0), cube.sum(axis=1)
    genes = list(adata.var_names)
    n = s1["NATMI"].reset_index(drop=True)
    share = lambda num, den: np.divide(num, den, out=np.zeros(np.shape(num)), where=np.asarray(den) != 0)
    rng = np.random.default_rng(0)
    B = N_PERMS
    for i in rng.choice(len(n), min(40, len(n)), replace=False):
        r = n.iloc[i]
        cs, ct = types.index(r.source), types.index(r.target)
        gl, gr = genes.index(r.ligand), genes.index(r.receptor)
        Ls, Rs = share(means[cs, gl], tot_obs[gl]), share(means[ct, gr], tot_obs[gr])
        s = r.magnitude
        assert s == pytest.approx(Ls * Rs, rel=1e-6, abs=1e-12)   # LIANA+'s means are float32
        Lp, Rp = share(cube[:, cs, gl], tot[:, gl]), share(cube[:, ct, gr], tot[:, gr])
        assert abs(r.p_null_a - np.mean(Lp * Rp >= s)) <= 1 / B + 1e-12     # a float32 tie may fall either side
        if s > 0:
            assert abs(r.p_L_null_b - (1 + np.sum(Lp * Rs >= s)) / (B + 1)) <= 1 / (B + 1) + 1e-12
            assert abs(r.p_R_null_b - (1 + np.sum(Ls * Rp >= s)) / (B + 1)) <= 1 / (B + 1) + 1e-12
        else:
            assert r.p_null_b == 1.0


def test_swap4_replaces_natmi_specificity(both):
    _, (lrs1, s1) = both
    V = lc.variants(lrs1, s1, nulls=NULLS)
    for k in NULLS:
        nm = _merge(V[f"swap4_{k}"], s1["NATMI"])
        np.testing.assert_array_equal(nm.natmi_p, nm[f"p_{k}"])
        cn = _merge(V[f"swap4_{k}"], s1["Connectome"])
        np.testing.assert_array_equal(cn.connectome_p, cn[f"p_{k}"])
    cp = _merge(V["swap4_null_a"], s1["CellPhoneDB"])
    np.testing.assert_array_equal(cp.cellphone_pvals, cp.p_null_a)        # null_a keeps the tools' own p
    cp = _merge(V["swap4_null_b"], s1["CellPhoneDB"])
    np.testing.assert_array_equal(cp.cellphone_pvals, cp.p_null_b)
