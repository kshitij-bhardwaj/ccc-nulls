"""Accuracy of the tail p-values (option 1 GPD vs option 2 Pearson III) against brute force.

    python -m tools.cellchat.b1_tail_accuracy --tissue mammary --n-lr 40 --B-ref 100000

Knijnenburg et al. (2009)'s validation design, on CellChat's own score and real Tabula Sapiens labels:
- Estimate: Null A with B = 1,000 (seed 1); every interaction with prob > 0 and fewer than 10
  permutations reaching it (where the tail models are used) gets a GPD and a Pearson III p.
- Reference: Null A again, independent permutations (seed 2), B_ref = 100,000 on the same
  interactions; p_ref = (1 + #) / (B_ref + 1), trusted where >= 10 permutations reach it (p_ref >= 1e-4).
- Interactions: n_lr ligand-receptor pairs drawn at random (seed 0) among those with such rows,
  all their sender-receiver rows (the cost of B_ref grows with the genes involved, not the rows).
- Report, per p_ref decade: median and IQR of log10(p_est / p_ref), and the share with
  p_est < p_ref / 2, i.e. anti-conservative by more than 2x (the error that inflates false calls).

Writes tools/cellchat/results/b1_tail_accuracy_<tissue>.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from tools.cellchat import adapter, nulls
from tools.cellchat.b1_calibration import log
from tools.cellchat.calibrate_null_a import TISSUES, load
from tools.cellchat.tail_pvalues import M_MIN, MIN_TAIL, N_EXC, permutation_p

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="mammary", choices=list(TISSUES))
    ap.add_argument("--n-lr", type=int, default=40)
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--B-ref", type=int, default=100_000)
    args = ap.parse_args()
    adata = load(args.tissue)
    df = adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(adata, "label"), nboot=2, seed=1))
    est = nulls.null_a(adata, df, "label", B=args.B, seed=1, rng="numpy", keep_top=N_EXC)
    k, tie, prob = (est[c].to_numpy() for c in ("n_exceed_null_a", "n_tie_null_a", "score_null_a"))
    in_tail = (prob > 0) & (k + tie < M_MIN) & (est["perm_npos_null_a"].to_numpy() >= MIN_TAIL)
    lrs = np.random.default_rng(0).permutation(np.unique(est.interaction_name[in_tail]))[:args.n_lr]
    sub = est.interaction_name.isin(lrs).to_numpy()
    log(f"{args.tissue}: {len(lrs)} L-R pairs, {sub.sum()} rows, {(sub & in_tail).sum()} in the tail")

    t0 = time.time()
    ref = nulls.null_a(adata, df[sub], "label", B=args.B_ref, seed=2, rng="numpy")
    log(f"reference B = {args.B_ref:,} in {time.time() - t0:.0f} s")
    k_ref = ref["n_exceed_null_a"].to_numpy()
    p_ref = (1 + k_ref) / (args.B_ref + 1)

    rows = in_tail[sub] & (k_ref >= M_MIN)                       # tail rows the reference resolves
    est_sub = est[sub].copy()
    est_sub.attrs = est.attrs
    p = {e: permutation_p(est_sub, e, seed=0)[0] for e in ("gpd", "pearson3")}
    p["permutation"] = ((k + tie + 1) / (args.B + 1))[sub]    # what Null A can say there: p <= ~1/100
    out = dict(dataset=TISSUES[args.tissue]["title"], B=args.B, B_ref=args.B_ref, n_lr=len(lrs),
               lrs=list(map(str, lrs)), n_rows=int(sub.sum()), n_compared=int(rows.sum()),
               n_tail_beyond_reference=int((in_tail[sub] & (k_ref < M_MIN)).sum()), by_decade={})
    for lo, hi in [(1e-4, 1e-3), (1e-3, 1e-2), (1e-4, 1e-2)]:
        r = rows & (p_ref >= lo) & (p_ref < hi)
        out["by_decade"][f"[{lo:g}, {hi:g})"] = {"n": int(r.sum())} | {
            e: dict(median_log10_err=float(np.median(np.log10(v[r] / p_ref[r]))),
                    iqr_log10_err=[float(x) for x in np.percentile(np.log10(v[r] / p_ref[r]), [25, 75])],
                    share_anticonservative_2x=float(np.mean(v[r] < p_ref[r] / 2)))
            for e, v in p.items() if r.any()}
    f = HERE / "results" / f"b1_tail_accuracy_{args.tissue}.json"
    f.write_text(json.dumps(out, indent=2))
    log("wrote", f)
    print(json.dumps(out["by_decade"], indent=1))


if __name__ == "__main__":
    main()
