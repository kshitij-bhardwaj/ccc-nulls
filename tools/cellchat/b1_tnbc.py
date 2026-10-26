"""Null B (B1) vs Null A on TNBC: the held-out test, run once with the settings frozen on Tabula Sapiens.

    python -m tools.cellchat.b1_tnbc

Design (each choice and why):
- Data: TNBC (Wu et al. 2021, Open Problems), all cells; truth = Open Problems' cytokine
  ligand -> target cell type pairs. Not used in any earlier B1 step.
- One CellChat pipeline run (its own candidate selection and score) and one Null A run
  (B = 1,000, keep_top for the tail): every procedure below post-processes the same output.
- Procedures (calls at alpha = 0.05):
    cellchat_p<0.05      CellChat as users run it: Null A p < 0.05, no correction.
    null_a_bh_scored     Null A, BH over CellChat's reported rows (prob > 0).
    null_a_bh            Null A, BH, on the same z / p and the same rows as B1 (theoretical null).
    b1                   B1, frozen settings (FROZEN, from b1_calibration on Tabula Sapiens).
    b1_<p-value>         B1 with the other tail model / without one (sensitivity).
- Ground truth: Open Problems' metrics with the procedure as a gate (a row keeps CellChat's prob
  only if called, as Open Problems gates CellPhoneDB by p); 95% CIs and the paired AUPRC difference
  from 1,000 resamples of the truth rows (match_null_a.bootstrap_metrics). Ranking by z alone is the
  same under both nulls (the empirical null is a monotone map of z), so the gate is where they differ.
- Disagreements: calls made by Null A only and by B1 only: how many, their prob and permutation
  z, and their precision against the truth among rows whose (ligand, target) the truth covers.
- Reproducibility: patients split in two halves (seeded); pipeline + Null A on each; Jaccard index
  of each procedure's call set (ligand-receptor, sender, receiver) between halves.

Writes tools/cellchat/results/b1_tnbc_wu2021.json.
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
from tools.cellchat.tail_pvalues import N_EXC, permutation_p

HERE = Path(__file__).resolve().parent
DATASET = "tnbc_wu2021"
ALPHA = 0.05
FROZEN = dict(pvalue="gpd", method="mle", pct0=0.25, family=None, fit_rows="nondegenerate")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def run(adata, groupby, B, seed):
    """Pipeline + Null A once; q-values of every procedure (call: q <= ALPHA)."""
    df = adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(adata, groupby), nboot=2, seed=1))
    res = nulls.null_a(adata, df, groupby, B=B, seed=seed, rng="numpy", keep_top=N_EXC)
    scored = res["score_null_a"].to_numpy() > 0
    q = {"cellchat_p<0.05": np.where(scored, res["p_null_a"].to_numpy(), 1.0)}
    qs = np.ones(len(res))
    qs[scored] = bh(res["p_null_a"].to_numpy()[scored])
    q["null_a_bh_scored"] = qs
    fits = {}
    for e in ("gpd", "pearson3", "permutation"):
        pv = permutation_p(res, e, seed)
        for name, method in [(f"null_a_bh_{e}", "theoretical"), (f"b1_{e}", FROZEN["method"])]:
            b = nulls.null_b(adata, df, groupby, null_a_result=res, seed=seed, pvalues=pv, **dict(FROZEN, pvalue=e, method=method))
            q[name] = b["q_null_b"].to_numpy()
            fits[name] = dict(b.attrs["null_b"]["fits"]["all"], how=b.attrs["null_b"]["how"])
            if name == f"b1_{FROZEN['pvalue']}":
                z = b["z_null_b"].to_numpy()
    q["null_a_bh"], q["b1"] = q.pop(f"null_a_bh_{FROZEN['pvalue']}"), q.pop(f"b1_{FROZEN['pvalue']}")
    return res, q, z, fits


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
    out = dict(dataset=DATASET, B=args.B, alpha=ALPHA, frozen=FROZEN, n_cells=int(adata.n_obs))

    res, q, z, fits = run(adata, groupby, args.B, seed=0)
    out["fits"] = fits
    called = {m: v <= ALPHA for m, v in q.items()}
    out["calls"] = {m: int(c.sum()) for m, c in called.items()}
    log("calls:", out["calls"])

    rows = res[["ligand_complex", "target"]]
    prob = res["score_null_a"].to_numpy()
    joined = {}
    for m in ["cellchat_prob_ungated", *q]:
        pred = opm.predictions(rows, prob, keys, p=None if m == "cellchat_prob_ungated" else q[m], alpha=ALPHA)
        gt = opm._join(pred, truth, keys)
        joined[m] = (gt.response.to_numpy().astype(int), gt.score.to_numpy())
    rnd = truth[list(keys)].copy()
    rnd["score"] = np.random.default_rng(1).random(len(rnd))           # Open Problems' random control
    gt = opm._join(rnd, truth, keys)
    joined["random"] = (gt.response.to_numpy().astype(int), gt.score.to_numpy())
    out["openproblems"] = bootstrap_metrics(joined, ref="b1", n_boot=args.metric_boot)
    for m, d in out["openproblems"]["methods"].items():
        log(f"  {m}: AUPRC {d['auprc']:.3f} {np.round(d['auprc_ci95'], 3).tolist()}, OR {d['odds_ratio']:.2f}",
            f"B1 - this {d['auprc_diff_ref_minus_this']:+.3f} {np.round(d['auprc_diff_ci95'], 3).tolist()}"
            if "auprc_diff_ci95" in d else "")

    # disagreements between Null A and B1 (same rows, same z)
    t = truth.set_index(list(keys)).response
    key = list(zip(res.ligand_complex.astype(str), res.target.astype(str)))
    resp = np.array([t.get(k, np.nan) for k in key], dtype=float)
    out["disagreement"] = {}
    for name, sel in [("both", called["null_a_bh"] & called["b1"]),
                      ("null_a_only", called["null_a_bh"] & ~called["b1"]),
                      ("b1_only", ~called["null_a_bh"] & called["b1"])]:
        covered = sel & np.isfinite(resp)
        out["disagreement"][name] = dict(
            n=int(sel.sum()), median_prob=float(np.median(prob[sel])) if sel.any() else None,
            median_z=float(np.median(z[sel])) if sel.any() else None,
            n_truth_covered=int(covered.sum()),
            truth_precision=float(resp[covered].mean()) if covered.any() else None)
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
            r, qh, _, _ = run(sub, groupby, args.B, seed=h + 1)
            sets.append({m: keys_of(r, v <= ALPHA) for m, v in qh.items()})
        out["reproducibility"] = dict(patients=[sorted(half), sorted(set(patients) - half)])
        for m in sets[0]:
            a, b = sets[0][m], sets[1][m]
            out["reproducibility"][m] = dict(n_half1=len(a), n_half2=len(b), n_both=len(a & b),
                                             jaccard=len(a & b) / max(len(a | b), 1))
        log("reproducibility:", {m: round(v["jaccard"], 3) for m, v in out["reproducibility"].items() if m != "patients"})

    f = HERE / "results" / f"b1_{DATASET}.json"
    f.write_text(json.dumps(out, indent=2, default=float))
    log("wrote", f)


if __name__ == "__main__":
    main()
