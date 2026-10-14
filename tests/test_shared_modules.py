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


def test_validate_null_accepts_same_rows_with_p():
    df = schema.validate(_rows())
    out = schema.validate_null(df, df.assign(p_null_a=0.5), "a")
    assert (out.p_null_a == 0.5).all()


@pytest.mark.parametrize("breaker, msg", [
    (lambda d: d, "missing column"),
    (lambda d: d.assign(p_null_b=0.0), "outside"),
    (lambda d: d.assign(p_null_b=np.nan), "missing or non-numeric"),
    (lambda d: d.iloc[::-1].assign(p_null_b=0.5), "rows differ"),
    (lambda d: d.iloc[1:].assign(p_null_b=0.5), "rows returned"),
])
def test_validate_null_rejects_bad_output(breaker, msg):
    df = schema.validate(_rows())
    with pytest.raises(ValueError, match=msg):
        schema.validate_null(df, breaker(df), "b")


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
