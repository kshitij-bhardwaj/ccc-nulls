"""Check that our CellChat Null A reproduces CellChat's permutation test on a dataset (Python only).

    python -m tools.cellchat.match_null_a --quick          # 8,000 TNBC cells
    python -m tools.cellchat.match_null_a                  # all 42,512 TNBC cells

1. R's generator: ``cellchat_py.RRandom`` must reproduce known R outputs, so the permutations
   are the ones CellChat draws with ``set.seed(seed.use)``.
2. CellChat, as ported line by line in ``cellchat_py`` (computeCommunProb, nboot = 100,
   seed.use = 1: its defaults), gives the reference prob and p-values.
3. ``nulls.null_a`` (the vectorised Null A) with the same seed must match them exactly: prob,
   p-values, BH q-values and significant calls.
4. ``nulls.null_a`` with independent numpy permutations (B = 100, 1,000) must agree up to
   Monte-Carlo error, i.e. the null distribution is the same.

Writes tools/cellchat/results/null_a_match_<dataset>_<tag>.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy import stats as sstats

from ccc_nulls import datasets
from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat.cellchat_py import RRandom

HERE = Path(__file__).resolve().parent
ALPHA = 0.05
# outputs of R >= 3.6 (default RNG kinds)
R_REFERENCE = {
    "set.seed(1); runif(3)": ([0.2655087, 0.3721239, 0.5728534], lambda: RRandom(1).runif(3).round(7).tolist()),
    "set.seed(123); runif(1)": ([0.2875775], lambda: RRandom(123).runif(1).round(7).tolist()),
    "set.seed(1); sample(10)": ([9, 4, 7, 1, 2, 5, 3, 10, 6, 8], lambda: RRandom(1).sample_int(10).tolist()),
    "set.seed(42); sample(10)": ([1, 5, 10, 8, 2, 4, 6, 9, 7, 3], lambda: RRandom(42).sample_int(10).tolist()),
    "set.seed(123); sample(10)": ([3, 10, 2, 8, 6, 9, 1, 7, 5, 4], lambda: RRandom(123).sample_int(10).tolist()),
}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def calls_table(ours, theirs):
    return dict(both=int((ours & theirs).sum()), only_ours=int((ours & ~theirs).sum()),
                only_cellchat=int((~ours & theirs).sum()), neither=int((~ours & ~theirs).sum()))


def compare(p_ref, p_ours):
    """p-values, BH q-values and calls of ours vs CellChat's (over the testable rows)."""
    qn, qo = bh(p_ref), bh(p_ours)
    return dict(
        n_rows=int(len(p_ref)),
        p_identical=int((p_ours == p_ref).sum()), p_max_abs_diff=float(np.abs(p_ours - p_ref).max()),
        p_mean_abs_diff=float(np.abs(p_ours - p_ref).mean()),
        p_spearman=float(sstats.spearmanr(p_ours, p_ref).statistic),
        q_identical=int((qo == qn).sum()), q_max_abs_diff=float(np.abs(qo - qn).max()),
        calls_p_lt_0_05=calls_table(p_ours < ALPHA, p_ref < ALPHA),
        calls_bh_q_lt_0_05=calls_table(qo < ALPHA, qn < ALPHA),
    )


def significance(p, B):
    q = bh(p)
    return dict(p_lt_0_05=int((p < ALPHA).sum()), p_eq_0=int((p == 0).sum()), bh_q_lt_0_05=int((q < ALPHA).sum()),
                bh_q_lt_0_10=int((q < 0.10).sum()), min_q=float(q.min()), p_resolution=1 / B)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="tnbc_wu2021")
    ap.add_argument("--quick", action="store_true", help="subsample 8,000 cells")
    ap.add_argument("--nboot", type=int, default=100, help="CellChat's nboot (its default: 100)")
    ap.add_argument("--seed", type=int, default=1, help="CellChat's seed.use (its default: 1)")
    ap.add_argument("--own-B", default="100,1000", help="independent-permutation runs of our Null A")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    tag = args.tag or ("quick" if args.quick else "full")
    groupby = datasets.DATASETS[args.dataset]["groupby"]
    out = dict(dataset=args.dataset, tag=tag, cellchat=adapter.CELLCHAT_REF, nboot=args.nboot, seed=args.seed)

    out["r_rng"] = {k: dict(expected=e, ours=f(), identical=f() == e) for k, (e, f) in R_REFERENCE.items()}
    log("R generator reproduces R:", all(v["identical"] for v in out["r_rng"].values()))

    log("loading", args.dataset)
    adata, _ = datasets.load_openproblems(args.dataset)
    if args.quick:
        rng = np.random.default_rng(0)
        adata = adata[np.sort(rng.choice(adata.n_obs, min(8000, adata.n_obs), replace=False))].copy()
        vc = adata.obs[groupby].value_counts()
        adata = adata[adata.obs[groupby].isin(vc[vc >= 5].index)].copy()
        adata.obs[groupby] = adata.obs[groupby].astype(str).astype("category")
    out.update(n_cells=int(adata.n_obs), n_groups=int(adata.obs[groupby].nunique()))
    log(f"{adata.n_obs:,} cells, {out['n_groups']} groups")

    t0 = time.time()
    inputs = adapter.prepare_inputs(adata, groupby, datasets.data_dir() / "cellchat" / f"{args.dataset}_{tag}")
    out["inputs"] = dict(n_genes_signaling=len(inputs["genes"]), n_features_sig=len(inputs["features_sig"]),
                         n_lr=int(len(inputs["LRsig"])), secs=round(time.time() - t0, 1))
    log(f"inputs: {len(inputs['genes'])} signaling genes, {len(inputs['features_sig'])} over-expressed, "
        f"{len(inputs['LRsig'])} L-R pairs ({time.time() - t0:.0f} s)")

    log(f"CellChat computeCommunProb (port), nboot = {args.nboot}, seed.use = {args.seed}")
    t0 = time.time()
    native = adapter.run_native(inputs, nboot=args.nboot, seed=args.seed, log=log)
    out["cellchat_secs"] = round(time.time() - t0, 1)
    df = adapter.native_to_schema(native)
    testable = df.score.to_numpy() > 0
    p_ref = df.p_native.to_numpy()[testable]
    out.update(n_rows=int(len(df)), n_rows_prob_gt_0=int(testable.sum()),
               cellchat_significance=significance(p_ref, args.nboot))
    log(f"  {len(df):,} rows, {int(testable.sum()):,} with prob > 0, {out['cellchat_secs']:.0f} s")

    log("our Null A, CellChat's permutations (rng = 'cellchat', same seed)")
    t0 = time.time()
    res = nulls.null_a(adata, df, groupby, B=args.nboot, seed=args.seed, rng="cellchat")
    prob_o, prob_n = res.score_null_a.to_numpy(), df.score.to_numpy()
    out["same_permutations"] = dict(
        secs=round(time.time() - t0, 1),
        prob_identical=int((prob_o == prob_n).sum()), prob_max_abs_diff=float(np.abs(prob_o - prob_n).max()),
        prob_zero_pattern_identical=bool(np.array_equal(prob_o == 0, prob_n == 0)),
        **compare(p_ref, res.p_null_a_cellchat.to_numpy()[testable]))
    r = out["same_permutations"]
    log(f"  prob identical {r['prob_identical']:,}/{len(df):,}; p identical {r['p_identical']:,}/{r['n_rows']:,}; "
        f"q identical {r['q_identical']:,}/{r['n_rows']:,}")

    out["independent_permutations"] = {}
    for B in [int(b) for b in args.own_B.split(",") if b]:
        log(f"our Null A, independent numpy permutations, B = {B}")
        t0 = time.time()
        ind = nulls.null_a(adata, df, groupby, B=B, seed=0, rng="numpy")
        po = ind.p_null_a_cellchat.to_numpy()[testable]
        cmp = compare(p_ref, po)
        # two independent Monte-Carlo estimates of the same p: sd = sqrt(p (1 - p) (1/B + 1/nboot))
        pbar = (p_ref + po) / 2
        sd = np.sqrt(pbar * (1 - pbar) * (1 / B + 1 / args.nboot))
        z = np.abs(po - p_ref) / np.where(sd > 0, sd, np.inf)
        cmp.update(secs=round(time.time() - t0, 1), significance=significance(po, B),
                   frac_within_2_mc_sd=float((z <= 2).mean()), frac_within_3_mc_sd=float((z <= 3).mean()))
        out["independent_permutations"][f"B{B}"] = cmp
        log(f"  spearman {cmp['p_spearman']:.4f}, mean |dp| {cmp['p_mean_abs_diff']:.4f}, "
            f"within 2 MC sd {cmp['frac_within_2_mc_sd']:.3f}")

    f = HERE / "results" / f"null_a_match_{args.dataset}_{tag}.json"
    f.parent.mkdir(exist_ok=True)
    f.write_text(json.dumps(out, indent=2))
    log("wrote", f)


if __name__ == "__main__":
    main()
