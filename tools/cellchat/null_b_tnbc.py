"""Null B vs Null A on TNBC: the held-out test, run once with the settings frozen on Tabula Sapiens.

    python -m tools.cellchat.null_b_tnbc

Design (each choice and why):
- Data: TNBC (Wu et al. 2021, Open Problems), all cells; truth = Open Problems' cytokine
  ligand -> target cell type pairs. Not used in any earlier Null B step.
- One CellChat pipeline run (its own candidate selection and score); Null A and Null B on the same
  permutations (B = 1,000, seed 0). Frozen setting: delta = 0 (null_b_calibration.py).
- Procedures (calls at alpha = 0.05):
    cellchat_p<0.05     CellChat as users run it: Null A p < 0.05, no correction.
    null_a              Null A, BH over the testable rows (nulls.testable).
    null_b              Null B = max(p_L, p_R), BH over the same rows.
    null_a|reported_rows  Null A, BH over CellChat's reported rows (prob > 0): what users do.
    null_b|delta=<d>    interval null (sensitivity).
- Ground truth: Open Problems' metrics with each procedure as a gate (a row keeps CellChat's prob
  only if called, as Open Problems gates CellPhoneDB by p); 95% CIs and the paired AUPRC difference
  Null B - other from 1,000 resamples of the truth rows (match_null_a.bootstrap_metrics).
- Disagreements: calls of Null A only / Null B only / both: how many, which side fails Null B,
  median prob, and precision against the truth among rows whose (ligand, target) the truth covers.
- Reproducibility: patients split into two halves (seeded); pipeline + both nulls on each half;
  Jaccard index of each procedure's call set (L-R pair, sender, receiver) between the halves.

Writes tools/cellchat/results/null_b_tnbc_wu2021.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from ccc_nulls import datasets
from ccc_nulls.evaluation import openproblems as opm
from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat.match_null_a import bootstrap_metrics
from tools.cellchat.null_b_calibration import DELTAS

HERE = Path(__file__).resolve().parent
DATASET = "tnbc_wu2021"
ALPHA = 0.05


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def run(adata, groupby, B, seed):
    """Pipeline, both nulls on the same permutations; q-values of every procedure (call: q <= ALPHA)."""
    df = adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(adata, groupby), nboot=2, seed=1))
    a = nulls.null_a(adata, df, groupby, B=B, seed=seed, rng="numpy")
    b = nulls.null_b(adata, df, groupby, B=B, seed=seed, rng="numpy", delta=DELTAS)
    scored = b.score_null_b.to_numpy() > 0
    family = nulls.testable(a, b)

    def bh_over(p, rows):
        q = np.ones(len(p))
        q[rows] = bh(p[rows])
        return q

    q = {"cellchat_p<0.05": np.where(scored, a.p_null_a.to_numpy(), 1.0),
         "null_a|reported_rows": bh_over(a.p_null_a.to_numpy(), scored),
         "null_a": bh_over(a.p_null_a.to_numpy(), family), "null_b": bh_over(b.p_null_b.to_numpy(), family),
         **{f"null_b|delta={d:g}": bh_over(b[f"p_null_b_delta={d:g}"].to_numpy(), family) for d in DELTAS[1:]}}
    return b, q


def keys_of(res, called):
    r = res[called]
    return set(zip(r.interaction_name.astype(str), r.source.astype(str), r.target.astype(str)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--metric-boot", type=int, default=1000)
    ap.add_argument("--no-halves", action="store_true")
    args = ap.parse_args()
    groupby = datasets.DATASETS[DATASET]["groupby"]
    sample_key = datasets.DATASETS[DATASET]["sample_key"]
    keys = datasets.DATASETS[DATASET]["truth_keys"]
    adata, truth = datasets.load_openproblems(DATASET)
    log(f"{DATASET}: {adata.n_obs:,} cells, {adata.obs[groupby].nunique()} groups")
    out = dict(dataset=DATASET, B=args.B, alpha=ALPHA, frozen=dict(delta=0.0), n_cells=int(adata.n_obs))

    res, q = run(adata, groupby, args.B, seed=0)
    called = {m: v <= ALPHA for m, v in q.items()}
    out["n_rows"], out["n_scored"] = int(len(res)), int((res.score_null_b.to_numpy() > 0).sum())
    out["calls"] = {m: int(c.sum()) for m, c in called.items()}
    log("calls:", out["calls"])

    rows = res[["ligand_complex", "target"]]
    prob = res.score_null_b.to_numpy()
    joined = {}
    for m in ["cellchat_prob_ungated", *q]:
        pred = opm.predictions(rows, prob, keys, p=None if m == "cellchat_prob_ungated" else q[m], alpha=ALPHA)
        gt = opm._join(pred, truth, keys)
        joined[m] = (gt.response.to_numpy().astype(int), gt.score.to_numpy())
    rnd = truth[list(keys)].copy()
    rnd["score"] = np.random.default_rng(1).random(len(rnd))           # Open Problems' random control
    gt = opm._join(rnd, truth, keys)
    joined["random"] = (gt.response.to_numpy().astype(int), gt.score.to_numpy())
    out["openproblems"] = bootstrap_metrics(joined, ref="null_b", n_boot=args.metric_boot)
    for m, d in out["openproblems"]["methods"].items():
        log(f"  {m}: AUPRC {d['auprc']:.3f} {np.round(d['auprc_ci95'], 3).tolist()}, OR {d['odds_ratio']:.2f}",
            f"Null B - this {d['auprc_diff_ref_minus_this']:+.3f} {np.round(d['auprc_diff_ci95'], 3).tolist()}"
            if "auprc_diff_ci95" in d else "")

    t = truth.set_index(list(keys)).response
    resp = np.array([t.get(k, np.nan) for k in zip(res.ligand_complex.astype(str), res.target.astype(str))], float)
    side = res.limiting_side_null_b.to_numpy()
    out["disagreement"] = {}
    for name, sel in [("both", called["null_a"] & called["null_b"]),
                      ("null_a_only", called["null_a"] & ~called["null_b"]),
                      ("null_b_only", ~called["null_a"] & called["null_b"])]:
        cov = sel & np.isfinite(resp)
        out["disagreement"][name] = dict(
            n=int(sel.sum()), sender_side_fails=int((sel & (side == "L")).sum()),
            receiver_side_fails=int((sel & (side == "R")).sum()),
            median_prob=float(np.median(prob[sel])) if sel.any() else None,
            n_truth_covered=int(cov.sum()), n_truth_positive=int(np.nansum(resp[cov])),
            truth_precision=float(resp[cov].mean()) if cov.any() else None)
    log("disagreement:", out["disagreement"])

    if not args.no_halves:
        patients = np.array(sorted(adata.obs[sample_key].astype(str).unique()))
        half = set(np.random.default_rng(0).permutation(patients)[: len(patients) // 2])
        in_half = adata.obs[sample_key].astype(str).isin(half).to_numpy()
        sets = []
        for h, mask in enumerate([in_half, ~in_half]):
            sub = adata[mask].copy()
            vc = sub.obs[groupby].value_counts()
            sub = sub[sub.obs[groupby].isin(vc[vc >= 5].index)].copy()
            sub.obs[groupby] = sub.obs[groupby].astype(str).astype("category")
            log(f"half {h + 1}: {sub.n_obs:,} cells")
            r, qh = run(sub, groupby, args.B, seed=h + 1)
            sets.append({m: keys_of(r, v <= ALPHA) for m, v in qh.items()})
        out["reproducibility"] = dict(patients=[sorted(half), sorted(set(patients) - half)])
        for m in sets[0]:
            a, b = sets[0][m], sets[1][m]
            out["reproducibility"][m] = dict(n_half1=len(a), n_half2=len(b), n_both=len(a & b),
                                             jaccard=len(a & b) / max(len(a | b), 1))
        log("reproducibility:", {m: round(v["jaccard"], 3) for m, v in out["reproducibility"].items() if m != "patients"})

    f = HERE / "results" / f"null_b_{DATASET}.json"
    f.write_text(json.dumps(out, indent=2, default=float))
    log("wrote", f)


if __name__ == "__main__":
    main()
