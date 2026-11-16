"""AMENDMENT_9's measures: the 1% share equals evaluate.top's false share, ties and all; the curve ends at the share
of all rows; auroc_false is 1 when every false row is ranked last."""
import numpy as np
import pandas as pd
import pytest

pytest.importorskip("liana")
from benchmarks.liana_final import amendment2 as A2  # noqa: E402
from benchmarks.liana_final import amendment9 as A9  # noqa: E402
from benchmarks.liana_final.evaluate import top  # noqa: E402


@pytest.fixture
def tbl():
    rng = np.random.default_rng(0)
    n = 1234
    d = {f"specificity_rank__{v}": rng.integers(0, 300, n) / 300 for v in ("liana6", "swap_null_b")}   # many ties
    return pd.DataFrame(d), rng.random(n) < 0.6


def test_one_percent_is_evaluate_top(tbl):
    t, f = tbl
    m = A9.cutoff_measures(t, f)
    for v in A2.variants_in(t):
        for c in A9.CUTS:
            assert m["false_share_at"][v][A9.key(c)] == pytest.approx(f[top(t, v, c)].mean(), abs=1e-15)
        assert m["false_share_curve"][v][-1] == pytest.approx(f.mean())


def test_auroc_false_direction(tbl):
    t, f = tbl
    t = t.assign(specificity_rank__liana6=np.where(f, 1.0, 0.0) + np.arange(len(t)) * 1e-9)
    assert A9.cutoff_measures(t, f)["auroc_false"]["liana6"] == pytest.approx(1.0)
