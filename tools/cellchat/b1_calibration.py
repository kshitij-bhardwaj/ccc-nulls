"""Null B (B1) calibration on Tabula Sapiens: the primary, statistical claim. TNBC is not touched.

    python -m tools.cellchat.b1_calibration --tissue all            # real labels + 20 full shuffles each
    python -m tools.cellchat.b1_calibration --tissue mammary --B-sensitivity 500,5000 --shuffles 3

Design (each choice and why):
- Known-null data: every cell's label shuffled once (seeds 1000 + r, as in calibrate_null_a.py).
  No interaction can be real, so every call is false: FDP = 1 if anything is called, else 0, and the
  FDR of a procedure is the mean FDP over shuffles = P(at least one call). A procedure controls FDR
  at alpha if that rate is <= alpha (Benjamini & Hochberg 1995).
- 20 shuffles per tissue, two tissues (Mammary, Liver), so the rate is estimated across independent
  null data sets, with a Clopper-Pearson 95% interval.
- Same candidates, observed scores and permutations for every procedure: CellChat's pipeline runs once
  per data set (its own candidate selection and score), Null A once (B = 1,000), and every procedure
  post-processes that one Null A output.
- Two factors, crossed:
    p-value   gpd (option 1), pearson3 (option 2), permutation (no tail model), naive (ablation:
              standardized score as z), see tail_pvalues.py;
    null      theoretical N(0, 1) = Null A with Benjamini-Hochberg on the same rows; empirical
              (Efron) = B1.
  So "B1 vs Null A" changes one thing, the null, with everything else held fixed.
- Reference: Null A as CellChat users run it, BH on p = (1 + #) / (B + 1) over its reported rows (prob > 0).
- alpha in {0.01, 0.05, 0.10}.
- p-value calibration (the input to both nulls): on shuffled data, the share of fitted rows with
  p <= t, t in {0.05, 0.01, 1e-3, 1e-4}; a valid p gives about t.
- Sensitivity of B1 (chosen here and frozen before TNBC): fit method (mle, cm) x central interval
  pct0 (0.10, 0.25, 0.33), family (all, per sender, per receiver), fitting rows (non-degenerate,
  CellChat's scored rows), B (--B-sensitivity).
- Real labels are run too, to show how many calls each procedure makes when signal is present.

Writes tools/cellchat/results/b1_calibration_<tissue><tag>.json.
"""
from __future__ import annotations

import argparse
import itertools
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats

from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat.calibrate_null_a import TISSUES, load
from tools.cellchat.tail_pvalues import N_EXC, permutation_p

HERE = Path(__file__).resolve().parent
ALPHAS = (0.01, 0.05, 0.10)
P_LEVELS = (0.05, 0.01, 1e-3, 1e-4)
ENGINES = ("gpd", "pearson3", "permutation", "naive")
DEFAULT = dict(pvalue="gpd", method="mle", pct0=0.25, family=None, fit_rows="nondegenerate")
VARIANTS = (
    [dict(DEFAULT, pvalue=e, method=m) for e, m in itertools.product(ENGINES, ["theoretical", "mle"])]
    + [dict(DEFAULT, method=m, pct0=p) for m, p in itertools.product(["mle", "cm"], [0.10, 0.25, 0.33])
       if (m, p) != ("mle", 0.25)]
    + [dict(DEFAULT, family=f) for f in ("source", "target")]
    + [dict(DEFAULT, fit_rows="scored"), dict(DEFAULT, pvalue="pearson3", fit_rows="scored")]
)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def vname(v):
    return f"{v['pvalue']}|{v['method']}|pct0={v['pct0']}|family={v['family']}|fit={v['fit_rows']}"


def one_dataset(adata, B, seed):
    """CellChat pipeline + Null A once; then every procedure from that one result."""
    inputs = adapter.prepare_inputs(adata, "label")
    df = adapter.native_to_schema(adapter.run_native(inputs, nboot=2, seed=1))
    t0 = time.time()
    res = nulls.null_a(adata, df, "label", B=B, seed=seed, rng="numpy", keep_top=N_EXC)
    scored = res["score_null_a"].to_numpy() > 0
    nondeg = res["n_tie_null_a"].to_numpy() < B
    out = dict(n_lr=int(df.interaction_name.nunique()), n_rows=int(len(df)), n_scored=int(scored.sum()),
               n_nondegenerate=int(nondeg.sum()), null_a_seconds=round(time.time() - t0, 1), calls={},
               p_calibration={}, how={})
    q = bh(res["p_null_a"].to_numpy()[scored]) if scored.any() else np.array([])
    out["calls"]["null_a|cellchat-convention|bh|scored"] = {str(a): int((q <= a).sum()) for a in ALPHAS}
    pvals = {e: permutation_p(res, e, seed) for e in ENGINES if e != "naive"}
    for e, (p, how) in pvals.items():
        pe = p[nondeg]
        out["p_calibration"][e] = {str(t): float((pe <= t).mean()) for t in P_LEVELS}
        out["how"][e] = {h: int(((how == h) & nondeg).sum()) for h in ("permutation", "tail", "capped")}
    z = res["z_null_a"].to_numpy()
    zn = z[np.isfinite(z)]
    out["p_calibration"]["naive"] = {str(t): float((sstats.norm.sf(zn) <= t).mean()) for t in P_LEVELS}
    for v in VARIANTS:
        try:
            b = nulls.null_b_empirical(adata, df, "label", null_a_result=res, seed=seed, pvalues=pvals.get(v["pvalue"]), **v)
        except ValueError as e:                                  # e.g. central matching not concave
            out["calls"][vname(v)] = {"error": str(e)}
            continue
        qb = b["q_null_b"].to_numpy()
        fit = b.attrs["null_b"]["fits"]
        out["calls"][vname(v)] = {str(a): int((qb <= a).sum()) for a in ALPHAS}
        out["calls"][vname(v)]["fit"] = fit.get("all", {k: fit[k] for k in list(fit)[:3]})
    return out


def clopper_pearson(k, n):
    lo = sstats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = sstats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return [float(lo), float(hi)]


def summarise(shuffled):
    """Per procedure and alpha: FDR = P(any call) over shuffles, with a 95% CI, and mean calls;
    per p-value method: mean share of p <= t."""
    out = {"fdr": {}, "p_calibration": {}}
    for name in shuffled[0]["calls"]:
        rows = [s["calls"][name] for s in shuffled if "error" not in s["calls"][name]]
        out["fdr"][name] = {}
        for a in ALPHAS:
            R = np.array([r[str(a)] for r in rows])
            k = int((R > 0).sum())
            out["fdr"][name][str(a)] = dict(fdr=k / len(R), fdr_ci95=clopper_pearson(k, len(R)),
                                            mean_calls=float(R.mean()), n_datasets=len(R))
    for e in shuffled[0]["p_calibration"]:
        out["p_calibration"][e] = {str(t): float(np.mean([s["p_calibration"][e][str(t)] for s in shuffled]))
                                   for t in P_LEVELS}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="all", choices=[*TISSUES, "all"])
    ap.add_argument("--shuffles", type=int, default=20)
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--B-sensitivity", default="", help="comma list of B to rerun on (subset)")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    for tissue in (list(TISSUES) if args.tissue == "all" else [args.tissue]):
        adata = load(tissue)
        log(f"{tissue}: {adata.n_obs:,} cells, {adata.obs.label.nunique()} cell types")
        out = dict(dataset=TISSUES[tissue]["title"], B=args.B, shuffles=args.shuffles, alphas=ALPHAS,
                   default=vname(DEFAULT), variants=[vname(v) for v in VARIANTS])
        f = HERE / "results" / f"b1_calibration_{tissue}{args.tag}.json"
        for B in [args.B] + [int(b) for b in args.B_sensitivity.split(",") if b]:
            key = "" if B == args.B else f"B{B}_"
            log(f"  real labels, B = {B}")
            out[key + "real_labels"] = one_dataset(adata, B, seed=0)
            shuffled = []
            labels = adata.obs.label.to_numpy()
            for r in range(args.shuffles):
                sh = adata.copy()
                sh.obs["label"] = pd.Categorical(np.random.default_rng(1000 + r).permutation(labels))
                shuffled.append(one_dataset(sh, B, seed=r + 1))
                d = shuffled[-1]["calls"]
                log(f"  shuffle {r + 1}/{args.shuffles}, B = {B}: calls at 0.05: "
                    + ", ".join(f"{e}: A {d[vname(dict(DEFAULT, pvalue=e, method='theoretical'))]['0.05']} "
                                f"B1 {d[vname(dict(DEFAULT, pvalue=e))].get('0.05')}" for e in ENGINES))
                out[key + "shuffled"] = shuffled
                out[key + "summary"] = summarise(shuffled)
                f.write_text(json.dumps(out, indent=2, default=float))     # partial results survive a stop
        log("wrote", f)


if __name__ == "__main__":
    main()
