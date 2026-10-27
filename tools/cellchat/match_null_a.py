"""Check that our CellChat Null A reproduces CellChat's permutation test on TNBC (Python only).

    python -m tools.cellchat.match_null_a --quick          # 8,000 TNBC cells
    python -m tools.cellchat.match_null_a                  # all 42,512 TNBC cells

1. R's generator: ``cellchat_py.RRandom`` must reproduce known R outputs, so the permutations
   are the ones CellChat draws with ``set.seed(seed.use)``.
2. For each CellChat version (2.1.2, the latest release, which is the default; 2.2.0.9001, GitHub
   main; 1.0.0, the 2021 paper's), the line-by-line port (``cellchat_py``, CellChat's defaults:
   nboot = 100, seed.use = 1) gives the reference prob and p-values, and the vectorised Null A
   (``nulls.null_a``) with the same seed must match them exactly: prob, p, BH q and calls.
3. Release 2.1.2: ``null_a`` with independent numpy permutations (B = 100, 1,000) must agree with
   CellChat up to Monte-Carlo error, i.e. the null distribution is the same.
4. How 2.2.0.9001 and 1.0.0 differ from 2.1.2.
5. Open Problems ligand-target metrics (TNBC cytokine truth) with 95% bootstrap CIs over the
   truth rows, paired across methods, including LIANA's CellChat re-implementation.

Writes tools/cellchat/results/null_a_match_<dataset>_<tag>.json.
"""
from __future__ import annotations

import argparse
import gc
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats
from sklearn.metrics import auc, precision_recall_curve

from ccc_nulls import datasets
from ccc_nulls.evaluation import openproblems as opm
from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat import cellchat_py as cc
from tools.cellchat.cellchat_py import RRandom

HERE = Path(__file__).resolve().parent
ALPHA = 0.05
RELEASE = cc.VERSION
VERSIONS = [RELEASE, "2.2.0.9001", "1.0.0"]
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


# ---------------------------------------------------------------------- metrics with bootstrap CIs
def _auprc(y, s):
    precision, recall, _ = precision_recall_curve(y, s, pos_label=1)
    return float(auc(recall, precision))


def _odds_ratio(y, s, top_prop=0.05):
    """Open Problems' odds ratio, with its tie order (pandas sort_values, descending): with many
    tied scores the top 5% is decided by that order, so it must be the same."""
    top = np.zeros(len(y), bool)
    top[pd.Series(s).sort_values(ascending=False).index.to_numpy()[: int(len(y) * top_prop)]] = True
    r = y == 1
    tp, fp, fn, tn = (top & r).sum(), (top & ~r).sum(), (~top & r).sum(), (~top & ~r).sum()
    if fp * fn == 0:
        return float("nan") if tp * tn == 0 else float("inf")
    return float(tp * tn / (fp * fn))


def bootstrap_metrics(joined: dict, ref: str, n_boot=1000, seed=0) -> dict:
    """``joined``: method -> (response, score) aligned on the truth rows. 95% percentile CIs from
    resampling truth rows with replacement (paired: the same resamples for every method), and
    the paired AUPRC difference ``ref`` minus each other method."""
    y = joined[ref][0]
    idx = np.random.default_rng(seed).integers(0, len(y), (n_boot, len(y)))
    boot = {m: np.array([_auprc(r[i], s[i]) for i in idx]) for m, (r, s) in joined.items()}
    orb = {m: np.array([_odds_ratio(r[i], s[i]) for i in idx]) for m, (r, s) in joined.items()}
    out = dict(n_truth_rows=int(len(y)), n_positive=int(y.sum()), n_boot=n_boot, reference=ref, methods={})
    for m, (r, s) in joined.items():
        fin = orb[m][np.isfinite(orb[m])]
        d = dict(n_truth_rows_scored=int((s > 0).sum()),
                 auprc=_auprc(r, s), auprc_ci95=[float(x) for x in np.percentile(boot[m], [2.5, 97.5])],
                 odds_ratio=_odds_ratio(r, s),
                 odds_ratio_ci95=[float(x) for x in np.percentile(fin, [2.5, 97.5])] if len(fin) else None,
                 odds_ratio_nonfinite_resamples=int((~np.isfinite(orb[m])).sum()))
        if m != ref:
            diff = boot[ref] - boot[m]
            d.update(auprc_diff_ref_minus_this=_auprc(*joined[ref]) - d["auprc"],
                     auprc_diff_ci95=[float(x) for x in np.percentile(diff, [2.5, 97.5])],
                     frac_resamples_ref_not_better=float((diff <= 0).mean()))
        out["methods"][m] = d
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="tnbc_wu2021")
    ap.add_argument("--quick", action="store_true", help="subsample 8,000 cells")
    ap.add_argument("--nboot", type=int, default=100, help="CellChat's nboot (its default: 100)")
    ap.add_argument("--seed", type=int, default=1, help="CellChat's seed.use (its default: 1)")
    ap.add_argument("--own-B", default="100,1000", help="independent-permutation runs of our Null A (release)")
    ap.add_argument("--metric-boot", type=int, default=1000)
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    tag = args.tag or ("quick" if args.quick else "full")
    groupby = datasets.DATASETS[args.dataset]["groupby"]
    out = dict(dataset=args.dataset, tag=tag, cellchat_release=adapter.CELLCHAT_REF,
               versions={v: cc.VERSIONS[v][0] for v in VERSIONS}, nboot=args.nboot, seed=args.seed)

    out["r_rng"] = {k: dict(expected=e, ours=f(), identical=f() == e) for k, (e, f) in R_REFERENCE.items()}
    log("R generator reproduces R:", all(v["identical"] for v in out["r_rng"].values()))

    log("loading", args.dataset)
    adata, truth = datasets.load_openproblems(args.dataset)
    keys = datasets.DATASETS[args.dataset]["truth_keys"]
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
        f"{len(inputs['LRsig'])} L-R pairs")

    official = {}

    def joined_for(rows, score, p=None, name=None):
        pred = opm.predictions(rows, score, keys, p=p)
        if name:
            official[name] = {k: float(v) for k, v in opm.metrics(pred, truth, keys).items()}
        gt = opm._join(pred, truth, keys)
        return gt.response.to_numpy().astype(int), gt.score.to_numpy()

    scores, pvals, joined = {}, {}, {}
    rows = None                                    # rows are in the same order in every version
    out["same_permutations"] = {}
    for v in VERSIONS:
        log(f"CellChat {v} ({cc.VERSIONS[v][0]}): port, nboot = {args.nboot}, seed.use = {args.seed}")
        t0 = time.time()
        df = adapter.native_to_schema(adapter.run_native(inputs, nboot=args.nboot, seed=args.seed, version=v))
        secs_port = time.time() - t0
        t0 = time.time()
        res = nulls.null_a(adata, df, groupby, B=args.nboot, seed=args.seed, rng="cellchat", version=v)
        testable = df.score.to_numpy() > 0
        p_ref = df.p_native.to_numpy()
        prob_o, prob_n = res.score_null_a.to_numpy(), df.score.to_numpy()
        out["same_permutations"][v] = dict(
            secs_port=round(secs_port, 1), secs_null_a=round(time.time() - t0, 1),
            n_rows_all=int(len(df)), n_rows_prob_gt_0=int(testable.sum()),
            prob_identical=int((prob_o == prob_n).sum()), prob_max_abs_diff=float(np.abs(prob_o - prob_n).max()),
            p_identical_all_rows=int((res.p_null_a_cellchat.to_numpy() == p_ref).sum()),
            **compare(p_ref[testable], res.p_null_a_cellchat.to_numpy()[testable]),
            cellchat_significance=significance(p_ref[testable], args.nboot))
        r = out["same_permutations"][v]
        log(f"  prob identical {r['prob_identical']:,}/{len(df):,}; p identical {r['p_identical']:,}/{r['n_rows']:,}; "
            f"q identical {r['q_identical']:,}/{r['n_rows']:,}")
        scores[v], pvals[v] = prob_n, p_ref
        if rows is None:
            rows = df[["ligand_complex", "target", "ligand_genes", "receptor_genes", "co_A_genes", "co_I_genes"]].copy()
        joined[f"cellchat_{v}"] = joined_for(rows, prob_n, name=f"cellchat_{v}")
        if v == RELEASE:
            joined[f"cellchat_{v}_gated_by_null_a"] = joined_for(rows, prob_n, p=res.p_null_a_cellchat.to_numpy(),
                                                                 name=f"cellchat_{v}_gated_by_null_a")
            out["independent_permutations"] = {}
            for B in [int(b) for b in args.own_B.split(",") if b]:
                log(f"  our Null A, independent numpy permutations, B = {B}")
                t0 = time.time()
                ind = nulls.null_a(adata, df, groupby, B=B, seed=0, rng="numpy", version=v)
                po, pr = ind.p_null_a_cellchat.to_numpy()[testable], p_ref[testable]
                cmp = compare(pr, po)
                # two independent Monte-Carlo estimates of the same p: sd = sqrt(p (1 - p) (1/B + 1/nboot))
                pbar = (pr + po) / 2
                sd = np.sqrt(pbar * (1 - pbar) * (1 / B + 1 / args.nboot))
                z = np.abs(po - pr) / np.where(sd > 0, sd, np.inf)
                cmp.update(secs=round(time.time() - t0, 1), significance=significance(po, B),
                           frac_within_2_mc_sd=float((z <= 2).mean()), frac_within_3_mc_sd=float((z <= 3).mean()))
                out["independent_permutations"][f"B{B}"] = cmp
                log(f"    spearman {cmp['p_spearman']:.4f}, mean |dp| {cmp['p_mean_abs_diff']:.4f}, "
                    f"within 2 MC sd {cmp['frac_within_2_mc_sd']:.3f}")
                del ind
        del df, res
        gc.collect()

    # how the other versions differ from the release
    s0, p0 = scores[RELEASE], pvals[RELEASE]
    has_cofactor = (rows.ligand_genes.str.contains(";") | rows.receptor_genes.str.contains(";")
                    | (rows.co_A_genes != "") | (rows.co_I_genes != "")).to_numpy()
    out["vs_release"] = {}
    for v in VERSIONS[1:]:
        s1, p1 = scores[v], pvals[v]
        both = (s0 > 0) & (s1 > 0)
        diff = s1 != s0
        out["vs_release"][v] = dict(
            prob_identical=int((~diff).sum()), n_rows=int(len(s0)),
            prob_differs_with_complex_or_coreceptor=int((diff & has_cofactor).sum()),
            prob_max_abs_diff=float(np.abs(s1 - s0).max()),
            prob_gt_0_release=int((s0 > 0).sum()), prob_gt_0_this=int((s1 > 0).sum()), prob_gt_0_both=int(both.sum()),
            p_identical_both=int((p1[both] == p0[both]).sum()), p_this_ge_release_both=int((p1[both] >= p0[both]).sum()),
            calls_p_lt_0_05=dict(both=int(((p1 < ALPHA) & (p0 < ALPHA)).sum()),
                                 only_this=int(((p1 < ALPHA) & (p0 >= ALPHA)).sum()),
                                 only_release=int(((p1 >= ALPHA) & (p0 < ALPHA)).sum())))
        log(f"  {v} vs {RELEASE}:", {k: x for k, x in out["vs_release"][v].items() if not isinstance(x, dict)},
            out["vs_release"][v]["calls_p_lt_0_05"])

    log("LIANA's CellChat re-implementation (lr_probs), for the ground-truth comparison")
    liana = adapter.liana_reference(adata, groupby, n_perms=None)
    joined["liana_cellchat"] = joined_for(liana, liana.score.to_numpy(), name="liana_cellchat")
    rnd = truth[list(keys)].copy()
    rnd["score"] = np.random.default_rng(1).random(len(rnd))           # Open Problems' random control
    gt = opm._join(rnd, truth, keys)
    joined["random"] = (gt.response.to_numpy().astype(int), gt.score.to_numpy())
    assert all(np.array_equal(r, joined["random"][0]) for r, _ in joined.values())
    log(f"bootstrap ({args.metric_boot} paired resamples of {len(truth)} truth rows)")
    out["openproblems"] = bootstrap_metrics(joined, ref=f"cellchat_{RELEASE}", n_boot=args.metric_boot)
    for m, o in official.items():             # point estimates must equal ccc_nulls.evaluation.openproblems
        d = out["openproblems"]["methods"][m]
        d["matches_openproblems_metrics"] = bool(np.isclose(d["auprc"], o["auprc"]) and
                                                 (np.isclose(d["odds_ratio"], o["odds_ratio"]) or
                                                  (np.isinf(d["odds_ratio"]) and np.isinf(o["odds_ratio"]))))
    for m, d in out["openproblems"]["methods"].items():
        log(f"  {m}: AUPRC {d['auprc']:.3f} {np.round(d['auprc_ci95'], 3).tolist()}, OR {d['odds_ratio']:.2f}",
            f"diff {np.round(d['auprc_diff_ci95'], 3).tolist()}" if "auprc_diff_ci95" in d else "")

    f = HERE / "results" / f"null_a_match_{args.dataset}_{tag}.json"
    f.parent.mkdir(exist_ok=True)
    f.write_text(json.dumps(out, indent=2))
    log("wrote", f)


if __name__ == "__main__":
    main()
