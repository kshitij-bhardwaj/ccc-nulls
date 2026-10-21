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
