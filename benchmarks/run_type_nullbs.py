"""Score the cell-type-ID permutation Null B candidates against Null A on the named simulations.

Reads the CellPhoneDB rows written by ``python -m tools.cellphonedb.run_sims``
(data/sim/cellphonedb/<file>.csv.gz, which hold the score at every cell-type pair) and the
matching simulations (data/sim/<file>.h5ad, for the planted truth). For every Null B module
present on this branch it computes p-values, runs Benjamini-Hochberg over the rows it can
test, and scores the calls exactly like Null A.

    python benchmarks/run_type_nullbs.py
Writes benchmarks/results/type_nullbs.json and prints one line per (file, null).
"""
import importlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ccc_nulls.stats import bh

ROOT = Path(__file__).resolve().parents[1]
NULLS = {"T1_type_exact": "ccc_nulls.nullb_type_exact",
         "T2_type_restd": "ccc_nulls.nullb_type_restd",
         "T3_type_spec": "ccc_nulls.nullb_type_spec"}


def available_nulls():
    out = {}
    for name, module in NULLS.items():
        try:
            out[name] = importlib.import_module(module).p_values
        except ImportError:
            pass
    return out


def score_null(res, p, sim, resource, alpha=0.05):
    """Metrics for one p-value vector, via the same evaluation as Null A."""
    from tools.cellphonedb.run_sims import evaluate

    ok = np.isfinite(p)
    q = np.full(len(p), np.nan)
    q[ok] = bh(p[ok])
    m = evaluate(res.assign(p_rand=p, q_rand=q, tested=ok), sim, resource, alpha)
    m["min_p"] = float(np.nanmin(p)) if ok.any() else None
    return m


def main():
    import anndata as ad

    from ccc_nulls.simulate import load_pairs, sim_from_anndata

    resource = load_pairs(ROOT / "data" / "resources" / "cellphonedb_pairs.csv")
    nulls = available_nulls()
    results = []
    for rows in sorted((ROOT / "data" / "sim" / "cellphonedb").glob("named_*.csv.gz")):
        name = rows.name.replace(".csv.gz", "")
        sim = sim_from_anndata(ad.read_h5ad(ROOT / "data" / "sim" / f"{name}.h5ad"))
        res = pd.read_csv(rows)
        runs = {"null_A": res.p_rand.where(res.tested).to_numpy(float)}
        runs.update({k: f(res) for k, f in nulls.items()})
        for null, p in runs.items():
            m = score_null(res, p, sim, resource)
            m.update(file=name, null=null)
            results.append(m)
            print(f"{name} {null:14s} tested {m['tested']:6d} calls {m['calls']:5d} TP {m['tp']:3d} "
                  f"FP {m['fp']:5d} partial {m['partial_calls']:4d} FDR {m['fdr']:.2f} power {m['power']} "
                  f"null p<.05 {m['fpr_p05']:.3f} min p {m['min_p']:.1e}", flush=True)
    out = ROOT / "benchmarks" / "results" / "type_nullbs.json"
    out.write_text(json.dumps(results, indent=1, default=float))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
