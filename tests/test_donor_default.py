"""rank_aggregate permutes within donor by default when the data have several donors."""
import warnings

import numpy as np
import pandas as pd
import pytest

li = pytest.importorskip("liana")
pytest.importorskip("anndata")

from ccc_nulls import liana_consensus as lc  # noqa: E402

from test_liana_consensus import N_PERMS, SEED, adata  # noqa: E402,F401


@pytest.fixture
def with_donors(adata):  # noqa: F811
    a = adata.copy()
    a.obs["donor_id"] = pd.Categorical(np.random.default_rng(1).choice(["d1", "d2", "d3"], a.n_obs))
    return a


def test_donor_column_found(adata, with_donors):  # noqa: F811
    assert lc.donor_strata(adata, "type") is None                       # no donor column: global
    assert lc.donor_strata(with_donors, "type") == "donor_id"
    assert lc.donor_strata(with_donors, "type", None) is None           # asked for global
    with pytest.raises(KeyError):
        lc.donor_strata(with_donors, "type", "batch")


def test_one_type_per_donor_falls_back_to_global(adata):  # noqa: F811
    a = adata.copy()
    a.obs["donor_id"] = a.obs["type"].astype(str)                      # every donor is one cell type
    with pytest.warns(UserWarning, match="global permutations"):
        assert lc.donor_strata(a, "type") is None


def test_rank_aggregate_default_is_within_donor(with_donors):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        auto, ta = lc.rank_aggregate(with_donors.copy(), "type", n_perms=N_PERMS, seed=SEED, inplace=False, return_tools=True)
        given, tg = lc.rank_aggregate(with_donors.copy(), "type", n_perms=N_PERMS, seed=SEED, inplace=False,
                                      return_tools=True, strata="donor_id")
        glob, tn = lc.rank_aggregate(with_donors.copy(), "type", n_perms=N_PERMS, seed=SEED, inplace=False,
                                     return_tools=True, strata=None)
    assert auto.attrs["strata"] == "donor_id" and glob.attrs["strata"] is None
    assert auto.equals(given)
    assert not np.array_equal(ta["CellChat"].p_null_b.to_numpy(), tn["CellChat"].p_null_b.to_numpy())
