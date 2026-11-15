from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from ccc_nulls import composite, schema
from tools.natmi import adapter

PAIRS = pd.DataFrame({"ligand": ["L1", "L1", "L2", "L3"], "receptor": ["R1", "R2", "R1", "R3"]})


def toy(seed: int = 0, n: int = 120):
    """Cells x genes, float32 ln(1 + x)-like values with zeros; four types; L3 and R3 never expressed in
    type 'd' (and R3 nowhere at all), so some scores are 0. Only what the adapter reads."""
    rng = np.random.default_rng(seed)
    genes = ["L1", "L2", "L3", "R1", "R2", "R3", "G1"]
    X = rng.gamma(1.0, 1.0, (n, len(genes))) * (rng.random((n, len(genes))) < 0.6)
    labels = np.array(list("abcd"))[rng.integers(0, 4, n)]
    X[labels == "d", 2] = 0
    X[:, 5] = 0
    obs = pd.DataFrame({"label": pd.Categorical(labels)}, index=[f"cell{i}" for i in range(n)])
    return SimpleNamespace(X=sparse.csr_matrix(np.log1p(X).astype(np.float32)), var_names=pd.Index(genes),
                           obs=obs, n_obs=n)


def natmi_reference(adata, pairs) -> pd.DataFrame:
    """NATMI's own arithmetic, transcribed from ExtractEdges.py (ClusterAnnotateEM, SplitIntoSinalProteins,
    GenSingleCell2CellEdge) in pandas: positive-score rows only, as NATMI writes them."""
    em = pd.DataFrame(adata.X.toarray().T.astype(float), index=adata.var_names)
    ann = adata.obs.label.astype(str).to_numpy()
    mean = pd.concat({c: em.loc[:, ann == c].mean(axis=1) for c in sorted(set(ann))}, axis=1)
    lig, rec = mean.loc[sorted(set(pairs.ligand))], mean.loc[sorted(set(pairs.receptor))]
    sl = lig.div(lig.sum(axis=1), axis=0).fillna(0.0)
    sr = rec.div(rec.sum(axis=1), axis=0).fillna(0.0)
    rows = []
    for p in pairs.itertuples():
        for s in mean.columns:
            for t in mean.columns:
                w = lig.loc[p.ligand, s] * rec.loc[p.receptor, t]
                if w > 0:
                    rows.append((p.ligand, p.receptor, s, t, sl.loc[p.ligand, s] * sr.loc[p.receptor, t], w))
    return pd.DataFrame(rows, columns=["ligand", "receptor", "source", "target", "score", "expr_prod"])


def test_score_is_natmis_arithmetic():
    a = toy()
    df, meta = adapter.run(a, "label", B=5, pairs=PAIRS)
    ref = natmi_reference(a, PAIRS)
    f = adapter.fidelity(df, ref)
    assert f["identical"], f
    assert f["only_port"] == 0 and f["only_native"] == 0 and f["rows_compared"] == len(ref) > 0
    assert len(df) == len(PAIRS) * 4 * 4 and meta["n_types"] == 4
    assert (df.loc[df.receptor == "R3", "score"] == 0).all()


def test_fidelity_catches_a_changed_score_and_a_missing_row():
    a = toy()
    df, _ = adapter.run(a, "label", B=2, pairs=PAIRS)
    ref = natmi_reference(a, PAIRS)
    bad = ref.copy()
    bad.loc[0, "score"] *= 1 + 1e-6
    assert not adapter.fidelity(df, bad)["identical"]
    f = adapter.fidelity(df, ref.iloc[1:])
    assert not f["identical"] and f["only_port"] == 1


def test_sides_split_exactly():
    """score(summarize(c), summarize(c)) is the score; the ligand side reads only the sender argument."""
    a = toy()
    h = adapter.side_score(a, "label", pairs=PAIRS)
    side, codes, K = h["side"], h["codes"], len(h["levels"])
    M = side.summarize(codes, K)
    s = side.score(M, M)
    M2 = side.summarize(np.random.default_rng(1).permutation(codes), K)
    sc = side.score
    lig_cols, rec_cols = sc.li, sc.ri
    assert len(np.intersect1d(lig_cols, rec_cols)) == 0
    assert not np.array_equal(M2, M)
    mixed = M.copy()
    mixed[:, rec_cols] = M2[:, rec_cols]              # receptor genes changed in the *sender* argument only
    assert np.array_equal(sc(mixed, M), s)
    mixed = M.copy()
    mixed[:, lig_cols] = M2[:, lig_cols]              # ligand genes changed in the *receiver* argument only
    assert np.array_equal(sc(M, mixed), s)


def test_null_a_uses_null_bs_permutations():
    """Null A's permutations are default_rng(seed).permutation(codes) per iteration, as composite's."""
    a = toy()
    h = adapter.side_score(a, "label", pairs=PAIRS)
    B, seed = 7, 3
    res = adapter.null_a(h, B=B, seed=seed)
    side, codes, K = h["side"], h["codes"], len(h["levels"])
    obs = side.summarize(codes, K)
    s_obs = side.score(obs, obs)
    gen = np.random.default_rng(seed)
    k = np.zeros(s_obs.shape, int)
    for _ in range(B):
        m = side.summarize(gen.permutation(codes), K)
        k += side.score(m, m) >= s_obs
    assert np.array_equal(res["n_exceed"], k)
    assert np.allclose(res["p"], np.where(s_obs == 0, 1.0, (1 + k) / (B + 1)))
    # Null B's engine draws the same sequence: with the receiver side observed and the sender side
    # permuted by the same draws, its sender counts are recomputed here
    gen = np.random.default_rng(seed)
    kL = np.zeros(s_obs.shape, int)
    for _ in range(B):
        kL += side.score(side.summarize(gen.permutation(codes), K), obs) >= s_obs
    b = adapter.null_b(h, B=B, seed=seed)
    assert np.allclose(b["p_L"], np.where(s_obs == 0, 1.0, (1 + kL) / (B + 1)))


def test_run_columns_and_family():
    a = toy()
    df, meta = adapter.run(a, "label", B=20, pairs=PAIRS)
    for c in ("p_null_a", "p_null_b", "p_L_null_b", "p_R_null_b", "limiting_side_null_b", "family", "expr_prod"):
        assert c in df
    assert (df.p_null_b >= df[["p_L_null_b", "p_R_null_b"]].max(axis=1) - 1e-15).all()
    zero_everywhere = df.receptor == "R3"
    assert not df.family[zero_everywhere].any()        # no labelling makes R3 positive
    assert (df.p_null_a[df.score == 0] == 1).all() and (df.p_null_b[df.score == 0] == 1).all()
    assert meta["n_family"] == int(df.family.sum())
    schema.validate(df)


def test_null_b_is_composite_on_natmis_score():
    a = toy()
    h = adapter.side_score(a, "label", pairs=PAIRS)
    direct = composite.composite(h["side"], h["codes"], len(h["levels"]), B=15, seed=2, ties="count")
    via = adapter.null_b(h, B=15, seed=2)
    assert np.array_equal(direct["p"], via["p"])


def test_side_genes_split_ligands_and_receptors():
    pairs = pd.DataFrame({"ligand": ["A", "B", "C"], "receptor": ["C", "D", "E"]})
    lig, rec = adapter.side_genes(pairs)
    assert lig == {"A", "B"} and rec == {"D", "E"}


def test_prepare_inputs_round_trips(tmp_path):
    a = toy()
    inp = adapter.prepare_inputs(a, "label", tmp_path, pairs=PAIRS)
    em = pd.read_csv(inp["em"], index_col=0)
    ann = pd.read_csv(inp["ann"], index_col=0)
    assert list(em.index) == sorted(adapter.database_genes(PAIRS))
    assert list(em.columns) == list(ann.index) and len(ann) == a.n_obs
    X = pd.DataFrame(a.X.toarray().T, index=a.var_names).loc[em.index]
    assert np.array_equal(em.to_numpy().astype(np.float32), X.to_numpy())


def test_liana_reference_schema(small_adata):
    pytest.importorskip("liana")
    df = adapter.liana_reference(small_adata, "label")
    assert len(df) > 0
    assert (df.score_name == adapter.SCORE_NAME).all()
    schema.validate(df)
    if adapter.HAS_NATIVE_NULL:
        assert df.p_native.notna().all()


def test_native_pipeline_schema(small_adata, tmp_path):
    try:
        inputs = adapter.prepare_inputs(small_adata, "label", tmp_path)
        native = adapter.run_native(inputs, tmp_path)
    except adapter.NativeNotInstalled as e:
        pytest.skip(str(e))
    df = adapter.native_to_schema(native)
    schema.validate(df, require_p_native=adapter.HAS_NATIVE_NULL)
    port, _ = adapter.run(small_adata, "label", B=2)
    f = adapter.fidelity(port, df)
    assert f["identical"], f
