import numpy as np
import pandas as pd

from ccc_nulls.nullb_type_restd import p_values


def table(n_inter=300, K=8, seed=0, effect=None):
    """Scores (L_a + R_b) / 2 with random cell-type effects: no cell-type pair is special."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_inter):
        L, R = rng.normal(1, 0.3, K), rng.normal(1, 0.3, K)
        if effect is not None and i == 0:
            L[0] += effect
            R[1] += effect
        for a in range(K):
            for b in range(K):
                rows.append((f"L{i}", f"R{i}", f"T{a}", f"T{b}", (L[a] + R[b]) / 2))
    df = pd.DataFrame(rows, columns=["ligand_complex", "receptor_complex", "source", "target", "score_obs"])
    return df.assign(tested=True)


def test_calibrated_without_special_pairs():
    p = p_values(table())
    assert 0.03 < (p < 0.05).mean() < 0.08


def test_planted_pair_is_extreme():
    res = table(n_inter=5, effect=3.0)
    p = p_values(res)
    hit = (res.ligand_complex == "L0") & (res.source == "T0") & (res.target == "T1")
    mine = (res.ligand_complex == "L0").to_numpy()
    assert p[hit.to_numpy()][0] < 1e-3 and p[hit.to_numpy()][0] == p[mine].min()
