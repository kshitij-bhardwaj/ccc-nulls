import os

import numpy as np
import pandas as pd
import pytest

from ccc_nulls import schema
from tools.cellphonedb import adapter, cpdb_null
from tools.cellphonedb.test_cpdb_null import toy


def cpdb_database():
    """Path to a downloaded CellPhoneDB database zip, from the CPDB_DB environment variable."""
    path = os.environ.get("CPDB_DB")
    if not path:
        pytest.skip("set CPDB_DB=/path/to/cellphonedb.zip to run the native CellPhoneDB tests")
    pytest.importorskip("cellphonedb")
    return path


def test_liana_reference_schema(small_adata):
    pytest.importorskip("liana")
    df = adapter.liana_reference(small_adata, "label", n_perms=10)
    assert len(df) > 0
    assert (df.score_name == adapter.SCORE_NAME).all()
    schema.validate(df)
    if adapter.HAS_NATIVE_NULL:
        assert df.p_native.notna().all()


def test_native_pipeline_schema(small_adata, tmp_path):
    db = cpdb_database()
    inputs = adapter.prepare_inputs(small_adata, "label", tmp_path)
    native = adapter.run_native(inputs, tmp_path, database_path=db)
    df = adapter.native_to_schema(native)
    if adapter.NATIVE_MODE == "method1":
        assert df.p_native.isna().all()            # Method 1 has no native p-values
    schema.validate(df, require_p_native=adapter.NATIVE_MODE == "method2")
    assert (df.ligand == df.ligand_complex).any()        # simple partners named by gene, not uniprot


def test_method1_with_nulls(small_adata, tmp_path):
    db = cpdb_database()
    res = adapter.method1_with_nulls(small_adata, "label", tmp_path, database_path=db, B=20)
    schema.validate(res)
    assert ((res.score - res.score_obs).abs() <= adapter.SCORE_TOLERANCE).all()
    assert res.tested.any()
    assert (res.ligand == res.ligand_complex).any()
    assert res.p_rand[res.tested].between(1 / 21, 1).all() and res.p_rand[~res.tested].isna().all()


def test_resource_pairs():
    pairs = adapter.resource_pairs(cpdb_database())
    assert list(pairs.columns) == ["ligand", "receptor"] and len(pairs) > 300
    assert not pairs.duplicated().any()
    assert ("CXCL12", "CXCR4") in set(map(tuple, pairs.to_numpy()))


def test_prepare_inputs_writes_h5ad(tmp_path):
    ad = pytest.importorskip("anndata")
    from scipy import sparse
    X = sparse.random(6, 4, density=0.5, format="csr", random_state=0)
    a = ad.AnnData(X=X, obs=pd.DataFrame({"ct": list("aabbcc")}, index=[f"AAAC-{i}" for i in range(6)]),
                   var=pd.DataFrame(index=["CXCL12", "CXCR4", "IL6", "IL6R"]))
    inputs = adapter.prepare_inputs(a, "ct", tmp_path)          # dashes in cell IDs are fine
    back = ad.read_h5ad(inputs["counts"])
    np.testing.assert_allclose(back.X.toarray(), X.toarray().astype(np.float32))
    meta = pd.read_csv(inputs["meta"], sep="\t")
    assert list(meta.columns) == ["Cell", "cell_type"] and list(meta.Cell) == list(a.obs_names)
    a.obs["ct"] = ["a|b"] * 6
    with pytest.raises(ValueError, match="pair separator"):
        adapter.prepare_inputs(a, "ct", tmp_path)


def fake_native(null, flip_ids=("I2",), undirected=("I3",)):
    """A Method 1 means table built from our own scores, rounded like CellPhoneDB's."""
    rows = []
    for i, iid in enumerate(null.interaction_ids):
        rec_a, rec_b = (True, False) if iid in flip_ids else (False, True)
        if iid in undirected:
            rec_a = rec_b = False
        rows.append({"id_cp_interaction": iid, "partner_a": f"simple:Q{i}A", "partner_b": f"simple:Q{i}B",
                     "gene_a": f"{iid}A", "gene_b": f"{iid}B", "receptor_a": rec_a, "receptor_b": rec_b})
    meta = pd.DataFrame(rows)
    means = pd.concat([meta, pd.DataFrame(np.round(null.score, 3), columns=null.pairs)], axis=1)
    means.loc[means.id_cp_interaction == "I2", "partner_a"] = "complex:cx"      # partner_a of I2 is the complex
    return {"means": means, "deconvoluted": None}


def test_combine_orientation_complexes_and_bh():
    data, _, _ = toy()
    null = cpdb_null.null_a(data, B=100, seed=3)
    with pytest.warns(UserWarning, match="both or neither"):
        res = adapter.combine(fake_native(null), null)
    assert "I3" not in set(res.native_interaction_id)                          # non-directional dropped
    simple = res[res.native_interaction_id != "I2"]
    assert (simple.ligand == simple.ligand_complex).all() and (simple.receptor == simple.receptor_complex).all()
    # I2 = (cx, p1) with partner A (cx) the receptor: ligand p1 sits in the pair's second type
    i2 = res[res.native_interaction_id == "I2"].set_index("cpdb_pair")
    row = i2.loc["T|B"]
    assert (row.source, row.target, row.ligand) == ("B", "T", "P1")
    assert row.receptor in {"P3", "P4"} and row.receptor_complex == "cx"
    k = data.pairs.index("T|B")
    assert row.score_obs == pytest.approx(null.score[1, k])
    # BH runs over tested rows only
    t = res.tested.to_numpy(bool)
    assert res.q_rand[~t].isna().all() and res.q_rand[t].notna().all()
    from ccc_nulls.stats import bh
    np.testing.assert_allclose(res.q_rand[t], bh(res.p_rand[t]))
    assert res.p_native.isna().all()


def test_combine_rejects_mismatched_scores():
    data, _, _ = toy()
    null = cpdb_null.null_a(data, B=10, seed=0)
    native = fake_native(null, undirected=())
    native["means"].loc[0, null.pairs[0]] += 0.01
    with pytest.raises(RuntimeError, match="beyond rounding"):
        adapter.combine(native, null)


def test_complex_named_like_a_gene_stays_distinct():
    """CellPhoneDB has e.g. complex "GALR1" and gene GALR1; both GAL -> GALR1 rows must survive."""
    data, _, _ = toy()
    null = cpdb_null.null_a(data, B=5, seed=0)
    native = fake_native(null, flip_ids=(), undirected=())
    m = native["means"]
    m.loc[m.id_cp_interaction == "I1", ["gene_a", "gene_b"]] = ["GAL", "GALR1"]
    m.loc[m.id_cp_interaction == "I4", ["gene_a", "partner_b"]] = ["GAL", "complex:GALR1"]
    df = adapter.native_to_schema(native)
    rc = set(df[df.ligand == "GAL"].receptor_complex)
    assert rc == {"GALR1", "complex:GALR1"}
    other = df[df.native_interaction_id == "I2"]                      # complex "cx": no clash, no prefix
    assert (other.ligand_complex == "cx").all()


def test_native_to_schema_without_null():
    data, _, _ = toy()
    null = cpdb_null.null_a(data, B=5, seed=0)
    with pytest.warns(UserWarning):
        df = adapter.native_to_schema(fake_native(null))
    i2 = df[df.native_interaction_id == "I2"]
    assert (i2.receptor == "cx").all() and (i2.ligand == "I2B").all()          # complex named, gene for simple
    assert (i2.ligand_complex == "I2B").all() and not df.ligand_complex.str.startswith("Q").any()
    assert len(df) == 3 * len(null.pairs)
