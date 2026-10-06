"""Calibration benchmark: does each null control false discoveries at the nominal level?

Regimes S0 (fully permuted), S1 (structured null) and S2 (planted signal, LFC 0.5–2.0);
see ccc_nulls.simulate. Writes benchmarks/results/simulation.json.

python benchmarks/simulation_study.py --seeds 5 --B 10000      # ~20 min on a laptop
python benchmarks/simulation_study.py --seeds 1 --B 1000 --quick
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

import ccc_nulls as cn
from ccc_nulls.simulate import candidates_from_pairs, label_truth, simulate

OUT = Path(__file__).parent / "results"


def run(seed, regime, B, score="mean", alpha=0.05, plant_lfc=1.0):
    sim = simulate(seed, structured=regime != "S0", planted=regime == "S2",
                   plant_lfc=plant_lfc, shuffle_labels=regime == "S0")
    cands = candidates_from_pairs(sim.X, sim.genes, sim.labels, sim.pairs)
    t0 = time.time()
    res = cn.significance_layer(sim.X, sim.genes, sim.labels, cands, score=score, B=B, alpha=alpha,
                                resource_genes=set(sim.pairs.values.ravel()), seed=seed)
    truth = label_truth(res, sim)
    is_true, is_null = truth == "true", truth == "null"
    mu = sim.X.mean(0)
    dec = pd.qcut(pd.Series(mu).rank(method="first"), 10, labels=False).to_numpy()
    lig_dec = dec[pd.Index(sim.genes).get_indexer(res.ligand)]

    def summ(p, q):
        call = q < alpha
        fp, tp = int((call & is_null).sum()), int((call & is_true).sum())
        return dict(calls=int(call.sum()), fp=fp, tp=tp, partial_calls=int((call & (truth == "partial")).sum()),
                    fdr=fp / max(fp + tp, 1), power=tp / is_true.sum() if is_true.any() else None,
                    fpr_p05=float((p[is_null] < 0.05).mean()),
                    calls_by_lig_decile=[int((call & is_null & (lig_dec == d)).sum()) for d in range(10)])

    pA, pB, qA, qB = (res[c].to_numpy() for c in ("p_rand", "p_strat", "q_rand", "q_strat"))
    flag = res.flag.to_numpy()
    return dict(seed=seed, regime=regime, score=score, plant_lfc=plant_lfc, n_tests=len(res),
                n_true=int(is_true.sum()), n_partial=int((truth == "partial").sum()), n_null=int(is_null.sum()),
                secs=round(time.time() - t0, 1), A=summ(pA, qA), B=summ(pB, qB),
                overlap={k: int((flag == k).sum()) for k in ("both", "A_only", "B_only")},
                overlap_null={k: int(((flag == k) & is_null).sum()) for k in ("both", "A_only", "B_only")},
                p_hist_A=np.histogram(pA[is_null], bins=20, range=(0, 1))[0].tolist(),
                p_hist_B=np.histogram(pB[is_null], bins=20, range=(0, 1))[0].tolist())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--B", type=int, default=10000)
    ap.add_argument("--quick", action="store_true", help="S1 and S2 at LFC 1.5 only")
    args = ap.parse_args()
    jobs = []
    if args.quick:
        jobs = [("S1", "mean", 1.0), ("S2", "mean", 1.5)]
    else:
        jobs = [(r, s, 1.0) for s in ("mean", "product") for r in ("S0", "S1")]
        jobs += [("S2", "mean", lfc) for lfc in (0.5, 1.0, 1.5, 2.0)]
    res = []
    for regime, score, lfc in jobs:
        for seed in range(args.seeds):
            r = run(seed, regime, args.B, score, plant_lfc=lfc)
            res.append(r)
            print(f"{score:7s} {regime} lfc={lfc} seed={seed} tests={r['n_tests']} | "
                  f"A calls={r['A']['calls']} FDR={r['A']['fdr']:.2f} | "
                  f"B calls={r['B']['calls']} FDR={r['B']['fdr']:.2f} power={r['B']['power']}", flush=True)
    OUT.mkdir(exist_ok=True)
    (OUT / ("simulation_quick.json" if args.quick else "simulation.json")).write_text(json.dumps(res, indent=1))
