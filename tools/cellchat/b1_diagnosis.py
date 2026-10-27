"""Why B1 (Efron's empirical null) cannot be used with CellChat: the evidence, on Tabula Sapiens.

    python -m tools.cellchat.b1_diagnosis --tissue all

Efron's empirical null estimates the null distribution from the central bulk of the z-values, which
assumes that bulk is null (the "zero assumption"; Efron 2004, sec. 4; 2010, ch. 6, where the
estimate is called reliable when the null proportion is about 0.9 or more). This script measures,
on the real labels, whether CellChat's z-values have such a bulk:

1. Depletion. Over the testable rows, the share with prob > 0 against the share random labels give
   (the mean of #{prob_b > 0} / B). Fewer non-zero rows than random labels give means most
   interactions are non-null, in the direction of absence.
2. The z-values (probit of the GPD-tail permutation p, tail_pvalues.py) by class: prob = 0 rows,
   reported rows (prob > 0), and reported rows compared only with their non-zero permutations
   (p_c = P(prob_b >= prob | prob_b > 0), exactly uniform under the null given the row is reported).
3. The empirical-null fit under each standard setting (MLE and central matching, central interval
   pct0 = 0.10, 0.25, 0.33), on all testable rows and on the reported rows' conditional z: the fitted
   null and the number of calls at q <= 0.05. If the null were identifiable, the settings would agree.

Writes tools/cellchat/results/b1_diagnosis_<tissue>.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats as sstats

from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat.b1_calibration import log
from tools.cellchat.calibrate_null_a import TISSUES, load
from tools.cellchat.empirical_null import empirical_null
from tools.cellchat.tail_pvalues import MIN_TAIL, N_EXC, permutation_p

HERE = Path(__file__).resolve().parent
SETTINGS = [("mle", 0.10), ("mle", 0.25), ("mle", 0.33), ("cm", 0.10), ("cm", 0.25), ("cm", 0.33)]
Q = (0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)


def fits(z):
    out = {"theoretical": int((empirical_null(z, "theoretical", with_lfdr=False)["q"] <= 0.05).sum())}
    for m, p0 in SETTINGS:
        try:
            r = empirical_null(z, m, p0, with_lfdr=False)
            f = r["fit"]
            out[f"{m}|pct0={p0}"] = dict(delta0=f.delta0, sigma0=f.sigma0, pi0=f.pi0, calls=int((r["q"] <= 0.05).sum()))
        except ValueError as e:
            out[f"{m}|pct0={p0}"] = {"error": str(e)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="all", choices=[*TISSUES, "all"])
    ap.add_argument("--B", type=int, default=1000)
    args = ap.parse_args()
    for tissue in (list(TISSUES) if args.tissue == "all" else [args.tissue]):
        adata = load(tissue)
        df = adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(adata, "label"), nboot=2, seed=1))
        res = nulls.null_a(adata, df, "label", B=args.B, seed=0, rng="numpy", keep_top=N_EXC)
        B = args.B
        prob, k, tie, npos = (res[c].to_numpy() for c in ("score_null_a", "n_exceed_null_a", "n_tie_null_a", "perm_npos_null_a"))
        p, how = permutation_p(res, "gpd", seed=0)
        z = sstats.norm.isf(p)
        testable = tie < B                         # not all B + 1 values equal
        scored = prob > 0
        # reported rows against their non-zero permutations: exact where >= 10 reach prob, GPD tail beyond
        graded = scored & (npos >= MIN_TAIL)
        pc = (k + np.random.default_rng(0).uniform(size=len(res)) * (tie + 1)) / (npos + 1)
        tail = graded & (how == "tail")
        pc[tail] = p[tail] * (B + 1) / (npos[tail] + 1)
        zc = sstats.norm.isf(np.clip(pc[graded], 1e-300, 1))
        out = dict(dataset=TISSUES[tissue]["title"], B=B, n_rows=int(len(res)), n_testable=int(testable.sum()),
                   n_reported=int(scored.sum()),
                   depletion=dict(share_nonzero_observed=float(scored[testable].mean()),
                                  share_nonzero_random_labels=float((npos[testable] / B).mean())),
                   reported_rows_pvalue_source={h: int(((how == h) & scored).sum()) for h in ("permutation", "tail", "capped")},
                   z_quantiles={name: {str(q): float(np.quantile(v, q)) for q in Q}
                                for name, v in [("testable", z[testable]), ("prob_0", z[testable & ~scored]),
                                                ("reported", z[scored]), ("reported_conditional", zc)]},
                   null_a_bh_calls_reported=int((bh(res.p_null_a.to_numpy()[scored]) <= 0.05).sum()),
                   fits_testable=fits(z[testable]), fits_reported_conditional=fits(zc))
        f = HERE / "results" / f"b1_diagnosis_{tissue}.json"
        f.write_text(json.dumps(out, indent=2))
        log(tissue, "depletion", out["depletion"])
        for name in ("fits_testable", "fits_reported_conditional"):
            log(f"  {name}:", {s: (v["calls"] if isinstance(v, dict) and "calls" in v else v) for s, v in out[name].items()})
        log("wrote", f)


if __name__ == "__main__":
    main()
