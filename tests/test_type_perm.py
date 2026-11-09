import numpy as np
import pandas as pd
from scipy import stats

from ccc_nulls import type_perm
from ccc_nulls.type_perm import loo_t_pvalues   # tested_groups via the module: pytest would collect a bare test* name


def test_loo_t_matches_hand_computation():
    s = np.array([1.0, 2.0, 4.0, 3.0, 10.0])
    p = loo_t_pvalues(s)
    others = s[:-1]
    t = (10 - others.mean()) / (others.std(ddof=1) * np.sqrt(1 + 1 / 4))
    assert np.isclose(p[-1], stats.t.sf(t, 3))
    assert np.isnan(loo_t_pvalues(np.array([1.0, 2.0, 3.0]))).all()     # fewer than 3 references
    flat = loo_t_pvalues(np.array([1.0, 1.0, 1.0, 1.0, 2.0]))     # last value's references have no spread
    assert flat[-1] == 0 and (flat[:-1] > 0.5).all()


def test_tested_groups_positions():
    res = pd.DataFrame({"k": ["a", "a", "b", "b"], "tested": [True, False, True, True]}, index=[10, 11, 12, 13])
    got = {tuple(g.k.unique()): rows.tolist() for rows, g in type_perm.tested_groups(res, ["k"])}
    assert got == {("a",): [0], ("b",): [2, 3]}
