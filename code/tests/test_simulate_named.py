import numpy as np
import pandas as pd
import pytest

import ccc_nulls as cn
from ccc_nulls.simulate import (candidates_from_pairs, label_truth, load_pairs, restrict_truth,
                                simulate_named)

# Stand-in resource with real symbols and LIANA's complex notation (``A_B``)
PAIRS = pd.DataFrame([(f"LIG{i}", f"REC{i}") for i in range(80)]
                     + [("CXCL12", "CXCR4"), ("CXCL12", "ACKR3"), ("IL6", "IL6R"), ("TNF", "TNFRSF1A"),
                        ("TGFB1", "TGFBR1_TGFBR2"), ("CCL5", "CCR5"), ("CCL5", "CCR5"), ("CD8A", "CD8A")],
                     columns=["source", "target"])


def test_load_pairs_cleans_resource():
    p = load_pairs(PAIRS)
    assert list(p.columns) == ["ligand", "receptor"]
    assert not p.receptor.str.contains("_").any()                  # complex dropped
    assert not p.duplicated().any() and (p.ligand != p.receptor).all()
    assert len(p) == 80 + 5
    assert len(load_pairs(PAIRS, simple_only=False)) == 80 + 6


def test_simulate_named_structure_and_truth():
    sim = simulate_named(PAIRS, seed=1, n_background=1500, sizes=(200, 150, 150, 100), n_planted=20)
    assert sim.X.shape == (600, len(sim.genes))
    assert {"CXCL12", "CXCR4", "IL6R"} <= set(sim.genes) and "BG00000" in set(sim.genes)
    assert len(sim.truth) == 20
    planted_genes = [g for t in sim.truth for g in t[:2]]
    assert len(planted_genes) == len(set(planted_genes))           # each planted gene used once
    resource = set(map(tuple, load_pairs(PAIRS).to_numpy()))
    assert all((l, r) in resource for l, r, _, _ in sim.truth)
    again = simulate_named(PAIRS, seed=1, n_background=1500, sizes=(200, 150, 150, 100), n_planted=20)
    assert again.truth == sim.truth and np.array_equal(again.X, sim.X)


def test_restrict_truth_to_tool_resource():
    sim = simulate_named(PAIRS, seed=2, n_background=500, sizes=(100, 100, 100), n_planted=10)
    some = pd.DataFrame(sorted({(l, r) for l, r, _, _ in sim.truth})[:4], columns=["ligand", "receptor"])
    assert len(restrict_truth(sim, some)) == 4
    assert restrict_truth(sim, some) <= sim.truth


def test_plant_pairs_limits_planting():
    allowed = pd.DataFrame({"ligand": [f"LIG{i}" for i in range(0, 80, 2)] + ["CXCL12", "IL6", "ZZZ"],
                            "receptor": [f"REC{i}" for i in range(0, 80, 2)] + ["CXCR4", "IL6R", "YYY"]})
    sim = simulate_named(PAIRS, seed=6, n_background=300, sizes=(60, 60, 60), n_planted=25, plant_pairs=allowed)
    ok = set(map(tuple, allowed.to_numpy()))
    assert len(sim.truth) == 25 and all((l, r) in ok for l, r, _, _ in sim.truth)
    with pytest.raises(ValueError, match="plant_pairs"):
        simulate_named(PAIRS, seed=6, n_background=300, sizes=(60, 60, 60), plant_pairs=allowed.iloc[[-1]])


def test_named_regimes_behave_like_simulate():
    """S1: Null A over-calls and Null B is calibrated. S2: Null B finds planted pairs."""
    kw = dict(n_background=2500, sizes=(300, 250, 200, 150), n_planted=15)
    s1 = simulate_named(PAIRS, seed=3, planted=False, **kw)
    c1 = candidates_from_pairs(s1.X, s1.genes, s1.labels, s1.pairs)
    r1 = cn.significance_layer(s1.X, s1.genes, s1.labels, c1, B=300, seed=3,
                               resource_genes=set(s1.pairs.values.ravel()))
    assert (r1.p_rand < 0.05).mean() > 0.15
    assert (r1.p_strat < 0.05).mean() < 0.10
    s2 = simulate_named(PAIRS, seed=4, plant_lfc=2.0, **kw)
    c2 = candidates_from_pairs(s2.X, s2.genes, s2.labels, s2.pairs)
    r2 = cn.significance_layer(s2.X, s2.genes, s2.labels, c2, B=300, seed=4,
                               resource_genes=set(s2.pairs.values.ravel()))
    truth = label_truth(r2, s2)
    assert (truth == "true").sum() >= 10
    assert np.median(r2.p_strat[truth == "true"]) < np.median(r2.p_strat[truth == "null"])


def test_anndata_round_trip(tmp_path):
    ad = pytest.importorskip("anndata")
    from ccc_nulls.simulate import sim_from_anndata, sim_to_anndata
    sim = simulate_named(PAIRS, seed=5, n_background=300, sizes=(60, 60, 60), n_planted=5)
    a = sim_to_anndata(sim, regime="S2", seed=5)
    a.write_h5ad(tmp_path / "s.h5ad")
    back = sim_from_anndata(ad.read_h5ad(tmp_path / "s.h5ad"))
    assert back.truth == sim.truth and back.lig_up == sim.lig_up and list(back.genes) == list(sim.genes)
    np.testing.assert_allclose(back.X, sim.X.astype(np.float32), rtol=1e-6)
