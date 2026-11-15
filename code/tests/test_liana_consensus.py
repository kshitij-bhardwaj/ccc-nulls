"""Our nulls inside LIANA+'s consensus (ccc_nulls.liana_consensus): LIANA+'s own outputs are untouched,
and Null B from LIANA+'s permutations is the composite null on the same permutations."""
import warnings

import numpy as np
import pandas as pd
import pytest

li = pytest.importorskip("liana")
an = pytest.importorskip("anndata")
from scipy import sparse  # noqa: E402

from ccc_nulls import composite  # noqa: E402
from ccc_nulls import liana_consensus as lc  # noqa: E402

N_PERMS, SEED = 40, 1337


@pytest.fixture(scope="module")
def adata():
    res = li.rs.select_resource("consensus")
    single = res[~res.ligand.str.contains("_") & ~res.receptor.str.contains("_")].drop_duplicates("ligand").head(25)
    complexes = res[res.receptor.str.contains("_")].head(3)
    genes = sorted({g for x in pd.concat([single, complexes])[["ligand", "receptor"]].to_numpy().ravel() for g in x.split("_")})
    rng = np.random.default_rng(0)
    K, n = 4, 600
    labels = rng.integers(0, K, n)
    lam = rng.gamma(2, 1, (K, len(genes)))
    counts = rng.poisson(lam[labels])
    X = np.log1p(counts / counts.sum(1, keepdims=True).clip(1) * 1e4).astype(np.float32)
    a = an.AnnData(sparse.csr_matrix(X), obs=pd.DataFrame({"type": pd.Categorical([f"t{k}" for k in labels])},
                                                           index=[f"c{i}" for i in range(n)]))
    a.var_names = genes
    return a


@pytest.fixture(scope="module")
def run(adata):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED)


def _merge(x, y):
    return x.drop_duplicates(lc.KEY).merge(y.drop_duplicates(lc.KEY), on=lc.KEY, suffixes=("_x", "_y"))


def test_five_method_consensus_is_liana_rank_aggregate(adata, run):
    lrs, store = run
    for d in store.values():
        d["p_b2"] = d["p_null_b"]          # B2 needs a decoy universe; not tested here
    V = lc.variants(lrs, store, nulls=("null_a", "null_b"))
    ref = li.mt.rank_aggregate(adata.copy(), "type", n_perms=N_PERMS, seed=SEED, use_raw=False, inplace=False)
    m = _merge(V["liana5"], ref)
    assert len(m) == len(ref) == len(V["liana5"])
    for c in ("magnitude_rank", "specificity_rank", "cellphone_pvals"):
        np.testing.assert_array_equal(m[f"{c}_x"], m[f"{c}_y"])


@pytest.mark.parametrize("tool,method", [("CellPhoneDB", "cellphonedb"), ("CellChat", "cellchat")])
def test_captured_tool_is_liana_standalone(adata, run, tool, method):
    _, store = run
    mag, spec = lc.TOOLS[tool]
    s = getattr(li.mt, method)(adata.copy(), "type", n_perms=N_PERMS, seed=SEED, use_raw=False, inplace=False)
    m = _merge(store[tool].rename(columns={"magnitude": mag, "p_null_a": spec}), s)
    assert len(m) == len(s)
    np.testing.assert_array_equal(m[f"{mag}_x"], m[f"{mag}_y"])
    np.testing.assert_array_equal(m[f"{spec}_x"], m[f"{spec}_y"])


def test_null_b_from_liana_permutations_is_the_composite_null(adata, run):
    """CellChat (trimeans: no dependence on cell order): the same p as ccc_nulls.composite run on
    LIANA+'s own label permutations, row for row."""
    _, store = run
    d = store["CellChat"]
    X = adata.X.copy()
    X /= np.float32(X.max())                                  # LIANA+'s CellChat scaling (mat_max), in float32
    X = X.toarray()
    types = list(adata.obs["type"].cat.categories)
    codes = adata.obs["type"].cat.codes.to_numpy()
    rng = np.random.default_rng(SEED)
    perms = np.stack([np.argsort(rng.permutation(adata.n_obs)) for _ in range(N_PERMS)], axis=1)
    gi = {g: i for i, g in enumerate(adata.var_names)}
    li_, ri_ = d.ligand.map(gi).to_numpy(), d.receptor.map(gi).to_numpy()
    si, ti = d.source.map({t: i for i, t in enumerate(types)}).to_numpy(), d.target.map({t: i for i, t in enumerate(types)}).to_numpy()

    def trimeans(c, K):
        out = np.zeros((K, X.shape[1]))
        for k in range(K):
            q = np.quantile(X[c == k], q=[0.25, 0.75], axis=0)      # LIANA+'s _trimean, step by step
            out[k] = (q[0] + 2 * np.median(X[c == k], axis=0) + q[1]) / 4
        return out

    side = composite.SideScore(summarize=trimeans, score=lambda S, R: lc._cellchat(S[si, li_], R[ti, ri_]))
    ref = composite.composite(side, codes, len(types), permutations=perms)
    np.testing.assert_allclose(d.score_check, ref["score"], rtol=1e-6)   # observed trimeans: a separate float32 path in LIANA+
    np.testing.assert_array_equal(d.p_L_null_b, ref["p_L"])
    np.testing.assert_array_equal(d.p_R_null_b, ref["p_R"])


def test_variants_swap_only_the_p_value(run):
    lrs, store = run
    for d in store.values():
        d["p_b2"] = d["p_null_b"]
    V = lc.variants(lrs, store, nulls=("null_a", "null_b"))
    assert set(V) == {"liana5", "liana6", "extra_null_a", "swap_null_b", "extra_null_b"}
    m = _merge(V["liana6"], V["swap_null_b"])
    np.testing.assert_array_equal(m.lr_means_x, m.lr_means_y)               # magnitude untouched
    cp = _merge(V["swap_null_b"], store["CellPhoneDB"])
    np.testing.assert_array_equal(cp.cellphone_pvals, cp.p_null_b)          # specificity is ours
    assert {"lr_means__null_b", "cellchat_pvals__null_b"} <= set(V["extra_null_b"].columns)


def test_ranking_helpers():
    from benchmarks.liana_final.evaluate import auroc, top
    t = pd.DataFrame({"specificity_rank__v": [0.5, 0.01, 0.2, 0.9, 0.05]})
    np.testing.assert_array_equal(top(t, "v", frac=0.4), [False, True, False, False, True])
    pos, neg = np.array([0, 1, 0, 0, 1], bool), np.array([1, 0, 1, 1, 0], bool)
    assert auroc(t, "v", pos, neg) == 1.0
    assert auroc(t, "v", neg, pos) == 0.0


def test_one_stratum_is_liana_permutation(adata, run):
    """A single stratum draws LIANA+'s own permutations: every p is unchanged."""
    _, store = run
    strata = pd.Series("all", index=adata.obs_names)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _, s = lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED, strata=strata, tools_only=True)
    for tool in lc.TOOLS:
        m = _merge(store[tool], s[tool])
        for c in ("p_null_a", "p_null_b"):
            np.testing.assert_array_equal(m[f"{c}_x"], m[f"{c}_y"])


def test_permutations_stay_within_strata(adata):
    """Strata = the cell types themselves: no permutation moves a cell to another type, so every
    permuted mean is the observed mean."""
    import liana.method.sc._liana_pipe as pipe
    a = adata.copy()
    a.obs["@label"] = a.obs["type"]
    obs = np.stack([np.asarray(a.X[(a.obs.type == t).to_numpy()].mean(0)).ravel() for t in a.obs.type.cat.categories])
    with lc.stratified_permutations(a.obs["type"].astype(str)):
        cube = pipe._get_means_perms(a, n_perms=5, seed=SEED, agg_fun=np.mean, norm_factor=None, n_jobs=1, verbose=False)
    np.testing.assert_allclose(cube, np.broadcast_to(obs, cube.shape), rtol=1e-5)
    from liana.method._pipe_utils._get_mean_perms import _get_means_perms
    assert pipe._get_means_perms is _get_means_perms               # LIANA+'s own again after the block


def test_rank_aggregate_public_call(adata):
    a = adata.copy()
    X0 = a.X.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res, tools = lc.rank_aggregate(a, "type", n_perms=N_PERMS, seed=SEED, return_tools=True)
    assert (a.X != X0).nnz == 0                                  # the user's data is untouched
    assert a.uns["liana_res"] is res and res.specificity_rank.is_monotonic_increasing
    cp = _merge(res, tools["CellPhoneDB"])
    np.testing.assert_array_equal(cp.cellphone_pvals, cp.p_null_b)


def test_rows_in_blocks_give_the_same_numbers(adata, run, monkeypatch):
    """lazy_permutations: rows scored 7 at a time give exactly what one block gives (and LIANA+'s own
    outputs are checked against li.mt.* above, through the same path)."""
    monkeypatch.setattr(lc, "ROWS_AT_ONCE", 7)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lrs, store = lc.run(adata.copy(), "type", n_perms=N_PERMS, seed=SEED)
    for tool, d in store.items():
        assert len(d) > 7
        pd.testing.assert_frame_equal(d, run[1][tool][d.columns])     # other tests add p_b2 to the fixture
    for k, t in lrs.items():
        pd.testing.assert_frame_equal(t, run[0][k])
