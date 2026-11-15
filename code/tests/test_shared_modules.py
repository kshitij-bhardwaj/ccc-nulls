import numpy as np
import pandas as pd
import pytest

from ccc_nulls import datasets, schema
from ccc_nulls.evaluation import openproblems as opm


def _rows(n=6):
    return pd.DataFrame({"ligand": [f"L{i}" for i in range(n)], "receptor": [f"R{i}" for i in range(n)],
                         "source": ["A"] * n, "target": ["B"] * n, "score": np.linspace(0, 1, n)})


def test_schema_fills_defaults_and_orders_columns():
    out = schema.validate(_rows())
    assert list(out.columns[:5]) == schema.REQUIRED
    assert (out.ligand_complex == out.ligand).all() and out.p_native.isna().all()


@pytest.mark.parametrize("breaker, msg", [
    (lambda d: d.drop(columns="score"), "missing required"),
    (lambda d: d.assign(score=["x"] * len(d)), "numeric"),
    (lambda d: pd.concat([d, d.iloc[:1]]), "duplicated"),
    (lambda d: d.assign(p_native=2.0), "outside"),
])
def test_schema_rejects_bad_tables(breaker, msg):
    with pytest.raises(ValueError, match=msg):
        schema.validate(breaker(_rows()))


def test_schema_require_p_native():
    with pytest.raises(ValueError, match="p_native"):
        schema.validate(_rows(), require_p_native=True)


def test_schema_require_p_native_rejects_missing_values():
    with pytest.raises(ValueError, match="missing values in 'p_native'"):
        schema.validate(_rows().assign(p_native=np.nan), require_p_native=True)
    schema.validate(_rows().assign(p_native=0.5), require_p_native=True)


def test_openproblems_metrics_perfect_and_gated():
    truth = pd.DataFrame({"ligand": list("abcd"), "target": ["t"] * 4, "response": [1, 1, 0, 0]})
    rows = pd.DataFrame({"ligand_complex": list("abcd"), "target": ["t"] * 4})
    perfect = opm.predictions(rows, [4, 3, 2, 1], ["ligand", "target"])
    assert opm.auprc(perfect, truth, ["ligand", "target"]) == pytest.approx(1.0)
    gated = opm.predictions(rows, [4, 3, 2, 1], ["ligand", "target"], p=[0.01, 0.5, 0.01, 0.5])
    assert gated.set_index("ligand").score.to_dict() == {"a": 4, "b": 0, "c": 2, "d": 0}


def test_dataset_registry():
    for name in ("tnbc_wu2021", "allen_brain_atlas"):
        meta = datasets.DATASETS[name]
        assert meta["url"].startswith("https://") and len(meta["md5"]) == 32
        assert datasets.path(name).name == meta["filename"]
    with pytest.raises(ValueError):
        datasets.path("tabula_sapiens")


@pytest.mark.parametrize("name", ["tabula_sapiens_mammary", "tabula_sapiens_liver"])
def test_tabula_sapiens_registry(name, tmp_path, monkeypatch):
    """The calibration tissues are registered like tnbc_wu2021 (no download needed)."""
    meta = datasets.DATASETS[name]
    assert meta["url"].startswith("https://datasets.cellxgene.cziscience.com/") and meta["url"].endswith(".h5ad")
    assert len(meta["md5"]) == 32 and int(meta["md5"], 16) >= 0
    assert isinstance(meta["size"], int) and meta["size"] > 0
    assert meta["collection"] == datasets.TABULA_SAPIENS == datasets.DATASETS["tabula_sapiens"]["url"]
    assert meta["groupby"] == "label" and meta["truth_keys"] is None and meta["title"].startswith("Tabula Sapiens - ")
    monkeypatch.setenv("CCC_NULLS_DATA", str(tmp_path))
    assert datasets.path(name) == tmp_path / meta["filename"] == tmp_path / f"{name}.h5ad"


def test_load_tabula_sapiens_reads_raw_counts(tmp_path, monkeypatch):
    """ln(1 + CP10k) of raw/X with library sizes over all genes, gene subset, small types dropped."""
    ad = pytest.importorskip("anndata")
    from scipy import sparse

    counts = np.array([[1, 0, 3, 0], [0, 2, 0, 2], [4, 0, 0, 1], [0, 0, 5, 5], [2, 2, 2, 2]], dtype=np.float32)
    obs = pd.DataFrame({"cell_type": pd.Categorical(["A", "A", "B", "B", "C"]),
                        "donor_id": pd.Categorical(["d1", "d1", "d2", "d2", "d2"])},
                       index=[f"cell{i}" for i in range(5)])
    var = pd.DataFrame({"feature_name": pd.Categorical(["G1", "G2", "G3", "G3"])}, index=[f"ENSG{i}" for i in range(4)])
    raw = ad.AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
    a = raw.copy()
    a.raw = raw
    monkeypatch.setenv("CCC_NULLS_DATA", str(tmp_path))
    a.write_h5ad(datasets.path("tabula_sapiens_mammary"))

    out = datasets.load_tabula_sapiens("tabula_sapiens_mammary", genes=["G1", "G3"], min_cells=2)
    assert list(out.obs_names) == ["cell0", "cell1", "cell2", "cell3"]          # type C (1 cell) dropped
    assert list(out.var_names) == ["G1", "G3"]                                   # duplicate G3 keeps the first
    assert list(out.obs.label.cat.categories) == ["A", "B"] and "donor_id" in out.obs
    lib = counts.sum(1, keepdims=True)                                           # all genes, not the subset
    expect = np.log1p(1e4 * counts / lib)[:4][:, [0, 2]].astype(np.float32)
    np.testing.assert_allclose(out.X.toarray(), expect, rtol=1e-6)
    with pytest.raises(ValueError):
        datasets.load_tabula_sapiens("tnbc_wu2021")
