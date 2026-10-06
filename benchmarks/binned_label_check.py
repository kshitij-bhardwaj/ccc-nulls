"""Does "permute cell labels within expression bins" differ from Null A?

It does not: each gene still sees a uniform label permutation, so its marginal null is
Null A's. The bins have to act on gene identity, which is what Null B does.
Also shows that fully permuted data (S0) cannot tell the nulls apart.
Writes benchmarks/results/strategy_check.json.

python benchmarks/binned_label_check.py --seeds 3 --B 2000
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import ccc_nulls as cn
from ccc_nulls.simulate import candidates_from_pairs, simulate

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", type=int, default=3)
ap.add_argument("--B", type=int, default=2000)
args = ap.parse_args()

out = []
for regime in ("S0", "S1"):
    for seed in range(args.seeds):
        sim = simulate(seed, structured=regime != "S0", planted=False, shuffle_labels=regime == "S0")
        cands = candidates_from_pairs(sim.X, sim.genes, sim.labels, sim.pairs)
        types = np.array(sorted(set(sim.labels)))
        lab = np.searchsorted(types, sim.labels)
        gi = pd.Index(sim.genes)
        li, ri = gi.get_indexer(cands.ligand), gi.get_indexer(cands.receptor)
        si, ti = np.searchsorted(types, cands.source), np.searchsorted(types, cands.target)
        st = cn.gene_strata(sim.X)
        pool = np.ones(len(gi), bool)
        pool[gi.get_indexer(sim.pairs.values.ravel())] = False
        _, pA = cn.null_a(sim.X, lab, len(types), li, ri, si, ti, "mean", args.B, seed)
        _, pL = cn.null_b_binned_labels(sim.X, lab, len(types), li, ri, si, ti, "mean", args.B, seed + 2, strata=st)
        _, pB = cn.null_b(sim.X, lab, len(types), li, ri, si, ti, "mean", seed + 1, strata=st, decoy_pool=pool)
        row = dict(regime=regime, seed=seed, n_tests=len(li))
        for name, p in (("A", pA), ("B_literal", pL), ("B", pB)):
            row[name] = dict(fpr_p05=float((p < 0.05).mean()), calls=int((cn.bh(p) < 0.05).sum()))
        row["corr_A_literal"] = float(np.corrcoef(np.log(pA), np.log(pL))[0, 1])
        out.append(row)
        print(row, flush=True)

res_dir = Path(__file__).parent / "results"
res_dir.mkdir(exist_ok=True)
(res_dir / "strategy_check.json").write_text(json.dumps(out, indent=1))
