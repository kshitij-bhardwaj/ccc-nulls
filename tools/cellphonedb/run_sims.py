"""Null A on CellPhoneDB's own statistic, scored against the frozen named simulations.

For every data/sim/named_*.h5ad: run CellPhoneDB Method 1 + Null A (method1_with_nulls),
then compare the q < alpha calls with the planted truth that CellPhoneDB's database can see.

    python -m tools.cellphonedb.run_sims --B 1000
    python -m tools.cellphonedb.run_sims --files data/sim/named_S2_seed0.h5ad --B 200

Needs ccc_nulls.simulate.sim_from_anndata / restrict_truth (branch core/named-simulation
until merged) and CPDB_DB or --database pointing at cellphonedb.zip.

Writes tools/cellphonedb/results/sims_null_a.json (metrics per file) and
data/sim/cellphonedb/<file>.csv.gz (all rows, for the later Null A vs Null B comparison).

Scoring (simple interactions that CellPhoneDB tested only. Complexes are counted separately):
  true     planted (ligand, receptor, source, target) visible to CellPhoneDB
  partial  shares a planted ligand-in-source or receptor-in-target, but is not planted
  null     everything else. A call here is a false positive.
FDR = null calls / (null + true calls). Partial calls are reported separately and counted in
neither, as in benchmarks/simulation_study.py.
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def evaluate(res: pd.DataFrame, sim, resource: pd.DataFrame, alpha: float = 0.05) -> dict:
    """Null A error and power on one simulation (``res`` from method1_with_nulls)."""
    from ccc_nulls.simulate import label_truth, restrict_truth

    visible = restrict_truth(sim, resource)
    sim_v = dataclasses.replace(sim, truth=visible)
    simple = (res.ligand == res.ligand_complex) & (res.receptor == res.receptor_complex)
    r = res.loc[simple & res.tested.astype(bool)].reset_index(drop=True)
    if r.empty:
        raise ValueError("no tested simple interactions to score; check ligand/ligand_complex columns")
    truth = label_truth(r, sim_v)
    is_true, is_null, is_part = truth == "true", truth == "null", truth == "partial"
    call = (r.q_rand < alpha).to_numpy()
    p = r.p_rand.to_numpy()
    fp, tp = int((call & is_null).sum()), int((call & is_true).sum())

    # abundance: ligand's overall mean expression decile across all genes in the simulation
    mu = sim.X.mean(0)
    dec = pd.qcut(pd.Series(mu).rank(method="first"), 10, labels=False).to_numpy()
    lig_dec = dec[pd.Index(sim.genes).get_indexer(r.ligand)]
    null_calls_by_decile = [int((call & is_null & (lig_dec == d)).sum()) for d in range(10)]
    null_tests_by_decile = [int((is_null & (lig_dec == d)).sum()) for d in range(10)]

    return dict(
        rows=len(res), tested=int(res.tested.sum()), tested_simple=len(r),
        tested_complex=int((~simple & res.tested.astype(bool)).sum()),
        planted=len(sim.truth), planted_visible=len(visible), true_tested=int(is_true.sum()),
        calls=int(call.sum()), tp=tp, fp=fp, partial_calls=int((call & is_part).sum()),
        fdr=fp / max(fp + tp, 1),
        power=tp / len(visible) if visible else None,
        fpr_p05=float((p[is_null] < 0.05).mean()) if is_null.any() else None,
        null_calls_by_lig_decile=null_calls_by_decile,
        null_call_rate_by_lig_decile=[round(c / t, 4) if t else None
                                      for c, t in zip(null_calls_by_decile, null_tests_by_decile)],
        p_hist_null=np.histogram(p[is_null], bins=20, range=(0, 1))[0].tolist(),
    )


def main():
    import anndata as ad

    from ccc_nulls.simulate import sim_from_anndata
    from tools.cellphonedb.adapter import method1_with_nulls, resource_pairs

    ap = argparse.ArgumentParser()
    ap.add_argument("--files", nargs="*", type=Path, default=None)
    ap.add_argument("--database", default=os.environ.get("CPDB_DB"))
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--threshold", type=float, default=0.1)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", type=Path, default=ROOT / "tools" / "cellphonedb" / "results" / "sims_null_a.json")
    args = ap.parse_args()
    if not args.database:
        ap.error("pass --database /path/to/cellphonedb.zip or set CPDB_DB")
    files = args.files or sorted((ROOT / "data" / "sim").glob("named_*.h5ad"))
    if not files:
        ap.error("no data/sim/named_*.h5ad files; run benchmarks/make_named_sims.py first")

    resource = resource_pairs(args.database)
    rows_dir = ROOT / "data" / "sim" / "cellphonedb"
    rows_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for f in files:
        adata = ad.read_h5ad(f)
        sim = sim_from_anndata(adata)
        params = dict(adata.uns.get("sim_params", {}))
        t0 = time.time()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            res = method1_with_nulls(adata, "label", rows_dir / f.stem, database_path=args.database, B=args.B,
                                     threshold=args.threshold, seed=args.seed, alpha=args.alpha)
        res.to_csv(rows_dir / f"{f.stem}.csv.gz", index=False)
        m = evaluate(res, sim, resource, args.alpha)
        m.update(file=f.name, regime=str(params.get("regime", "?")), sim_seed=params.get("seed"),
                 B=args.B, threshold=args.threshold, alpha=args.alpha, secs=round(time.time() - t0, 1))
        results.append(m)
        print(f"{f.name}: tested {m['tested_simple']} simple (+{m['tested_complex']} complex) | "
              f"calls {m['calls']} (TP {m['tp']}, FP {m['fp']}, partial {m['partial_calls']}) | "
              f"FDR {m['fdr']:.2f} power {m['power'] if m['power'] is None else round(m['power'], 2)} | "
              f"null p<.05 {m['fpr_p05']} | {m['secs']}s", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
