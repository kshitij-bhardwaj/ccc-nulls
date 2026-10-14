"""Null-model comparison on the Open Problems / LIANA TNBC atlas (Wu et al. 2021).

Data:  openproblems_v1/tnbc_wu2021/log_cp10k (Open Problems S3), 42,512 cells, 29 types,
       10 patients, cytokine ground truth in uns["ccc_target"].
Tools: 7 liana-py methods, each tested on the statistic Open Problems ranks it by.
Nulls: N1 global label permutation, N2 patient-stratified permutation, N3 abundance-matched
       decoys, N4 unmatched decoys, N5 patient-replicate t-test, N6 LIANA's native p.
Evals: E1 Open Problems ligand-target AUPRC / odds ratio, E2 fake LR resource,
       E4 abundance profile, E5 shuffled-label sanity.

python benchmarks/tnbc_openproblems.py --B 1000 --B5 100           # full run
python benchmarks/tnbc_openproblems.py --B 50 --B5 10 --quick      # smoke test
"""
import argparse
import json
import pickle
import time
import warnings
from pathlib import Path

import liana as li
import numpy as np
import pandas as pd

import ccc_nulls as cn
from ccc_nulls import datasets
from ccc_nulls.evaluation import openproblems as opm
from ccc_nulls.multitool import (TOOLS, Candidates, StatsEngine, decoy_null, label_null,
                                 patient_replicate_null, tool_scores)

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "benchmarks" / "results" / "tnbc"
EXPR_PROP, ALPHA = 0.10, 0.05
NATIVE_P = {"CellPhoneDB": "cellphone_pvals", "GeometricMean": "gmean_pvals", "CellChat": "cellchat_pvals"}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ------------------------------------------------------------------ Open Problems metrics
KEYS = ["ligand", "target"]                      # ligand-target task (cytokine truth)


def op_metrics(pred: pd.DataFrame, truth: pd.DataFrame):
    return opm.metrics(pred, truth, KEYS)


def op_prediction(df: pd.DataFrame, score: np.ndarray, p=None):
    """Open Problems protocol: keep the tool score if p <= 0.05 (else 0), max per (ligand, target)."""
    return opm.predictions(df, score, KEYS, p=p, alpha=ALPHA)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1000, help="permutations / decoy draws for N1-N4")
    ap.add_argument("--B5", type=int, default=100, help="within-patient permutations for N5")
    ap.add_argument("--quick", action="store_true", help="subsample cells for a smoke test")
    ap.add_argument("--tools", default=",".join(TOOLS))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--nulls", default="N1,N2,N3,N4,N5")
    ap.add_argument("--no-e5", action="store_true")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    tools = args.tools.split(",")
    rng = np.random.default_rng(args.seed)
    OUT.mkdir(parents=True, exist_ok=True)
    tag = args.tag or ("quick" if args.quick else "full")
    liana_tag = "quick" if args.quick else "full"

    log("loading", datasets.path("tnbc_wu2021").name)
    adata, truth = datasets.load_openproblems("tnbc_wu2021")   # X = Open Problems' log-CP10k layer
    if args.quick:
        keep = rng.choice(adata.n_obs, 8000, replace=False)
        adata = adata[np.sort(keep)].copy()
        # LIANA drops cell types with < 5 cells; drop them up front so every null sees the same types
        vc = adata.obs["label"].value_counts()
        adata = adata[adata.obs["label"].isin(vc[vc >= 5].index)].copy()
        adata.obs["label"] = adata.obs["label"].astype(str).astype("category")
    patients = adata.obs["orig.ident"].astype(str).to_numpy()
    types = list(adata.obs.label.cat.categories)
    K = len(types)
    log(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes, {K} cell types, {len(set(patients))} patients, "
        f"truth {len(truth):,} rows ({truth.response.sum()} positive)")

    # ---------------- LIANA runs (cached)
    cache = OUT / f"liana_{liana_tag}.pkl"
    if cache.exists():
        liana_res = pickle.loads(cache.read_bytes())
    else:
        liana_res = {}
        for t in tools:
            m, col, _ = TOOLS[t]
            log("liana", m)
            kw = dict(groupby="label", resource_name="consensus", expr_prop=EXPR_PROP, use_raw=False,
                      verbose=False, key_added=m, seed=args.seed)
            if t in NATIVE_P:
                # LIANA's native permutations are slow (dense trimeans for CellChat, scipy gmean);
                # CellPhoneDB gets 1,000 for the main native-vs-N1 cross-check, the others 100.
                kw["n_perms"] = 1000 if (t == "CellPhoneDB" and not args.quick) else 100
            getattr(li.mt, m)(adata, **kw)
            liana_res[t] = adata.uns[m].copy()
        cache.write_bytes(pickle.dumps(liana_res))
    for t in tools:
        log(f"  {t}: {len(liana_res[t]):,} interactions")

    # ---------------- engines and validation against LIANA
    X = adata.X
    # LIANA drops all-zero genes before taking the global mean (prep_check_adata)
    nonempty = int((np.asarray(X.sum(0)).ravel() > 0).sum())
    mat_mean = float(X.sum(dtype="float64") / (X.shape[0] * nonempty))
    mat_max = float(X.max())
    genes = pd.Index(adata.var_names.astype(str))
    resource = li.resource.select_resource("consensus")
    res_genes = set(resource.ligand.str.split("_").explode()) | set(resource.receptor.str.split("_").explode())
    lr_idx = np.where(genes.isin(list(res_genes)))[0]
    lr_pos = {g: i for i, g in enumerate(genes[lr_idx])}
    all_pos = {g: i for i, g in enumerate(genes)}
    type_pos = {t: i for i, t in enumerate(types)}
    codes = adata.obs.label.cat.codes.to_numpy()

    log("stats engines")
    need = set().union(*(TOOLS[t][2] for t in tools)) | {"props"}
    S_obs = StatsEngine(X, lr_idx, mat_mean, mat_max).compute(codes, K, need)
    # whole-transcriptome statistics (for decoys) under real and fully shuffled labels, then free
    # the engine: three copies of a 42k x 28k matrix do not fit next to everything else in 8 GB
    codes_sh = np.random.default_rng(args.seed + 7).permutation(codes)
    eng_all = StatsEngine(X, np.arange(len(genes)), mat_mean, mat_max)
    S_all = eng_all.compute(codes, K, need)
    S_all_sh = eng_all.compute(codes_sh, K, need)
    del eng_all
    import gc; gc.collect()
    log("  done")

    cands = {t: Candidates(liana_res[t], lr_pos, type_pos) for t in tools}
    cands_all = {t: Candidates(liana_res[t], all_pos, type_pos) for t in tools}
    obs, validation = {}, {}
    for t in tools:
        c = cands[t]
        obs[t] = tool_scores(t, S_obs, c.li, c.ri, c.si, c.ti, mat_mean)
        ref = liana_res[t][TOOLS[t][1]].to_numpy(dtype=float)
        validation[t] = dict(max_abs_diff=float(np.nanmax(np.abs(obs[t] - ref))),
                             corr=float(np.corrcoef(obs[t], ref)[0, 1]))
        log(f"  validate {t}: max|ours - liana| = {validation[t]['max_abs_diff']:.2e}, r = {validation[t]['corr']:.6f}")

    # ---------------- fake LR resource (E2): matched pairs from a held-out half of background genes
    strata = cn.gene_strata(X)
    bg = np.setdiff1d(np.arange(len(genes)), [all_pos[g] for g in res_genes if g in all_pos])
    bg = bg[np.asarray((X[:, bg] > 0).sum(0)).ravel() > 0]
    half = rng.permutation(bg)
    decoy_pool = np.zeros(len(genes), bool); decoy_pool[half[: len(half) // 2]] = True
    fake_pool = np.zeros(len(genes), bool); fake_pool[half[len(half) // 2:]] = True
    props_all = S_all["props"]
    fake = {}
    for t in tools:
        df = liana_res[t]
        pairs = df[["ligand", "receptor"]].drop_duplicates()
        fl, fr = [], []
        for L, R in zip(pairs.ligand, pairs.receptor):
            pl = np.where(fake_pool & (strata == strata[all_pos[L]]))[0]
            pr = np.where(fake_pool & (strata == strata[all_pos[R]]))[0]
            fl.append(genes[rng.choice(pl)]); fr.append(genes[rng.choice(pr)])
        fp = pd.DataFrame({"ligand": fl, "receptor": fr}).drop_duplicates()
        li_, ri_ = fp.ligand.map(all_pos).to_numpy(), fp.receptor.map(all_pos).to_numpy()
        rows = []
        for a, b, L, R in zip(li_, ri_, fp.ligand, fp.receptor):
            xs = np.where(props_all[:, a] >= EXPR_PROP)[0]
            ys = np.where(props_all[:, b] >= EXPR_PROP)[0]
            rows += [(L, R, types[x], types[y]) for x in xs for y in ys]
        fdf = pd.DataFrame(rows, columns=["ligand", "receptor", "source", "target"])
        fdf["ligand_complex"], fdf["receptor_complex"] = fdf.ligand, fdf.receptor
        fake[t] = fdf
        log(f"  fake resource {t}: {len(fp)} pairs -> {len(fdf):,} interactions")

    # fake rows need their own engine columns (they are background genes)
    fake_genes = sorted(set(pd.concat([fake[t].ligand for t in tools])) | set(pd.concat([fake[t].receptor for t in tools])))
    comb_idx = np.concatenate([lr_idx, [all_pos[g] for g in fake_genes]])
    comb_pos = {g: i for i, g in enumerate(genes[comb_idx])}
    eng_comb = StatsEngine(X, comb_idx, mat_mean, mat_max)
    S_comb = eng_comb.compute(codes, K, need)
    both = {t: pd.concat([liana_res[t].assign(is_fake=False), fake[t].assign(is_fake=True)], ignore_index=True) for t in tools}
    cands_b = {t: Candidates(both[t], comb_pos, type_pos) for t in tools}
    cands_b_all = {t: Candidates(both[t], all_pos, type_pos) for t in tools}
    obs_b = {t: tool_scores(t, S_comb, c.li, c.ri, c.si, c.ti, mat_mean) for t, c in cands_b.items()}

    # ---------------- nulls on real + fake rows
    P = {}
    nulls = args.nulls.split(",")
    if "N1" in nulls:
        log(f"N1 global permutation, B={args.B}")
        P["N1"] = label_null(eng_comb, codes, K, cands_b, obs_b, args.B, rng, None, tools, log_every=max(args.B // 5, 1))
    if "N2" in nulls:
        log(f"N2 patient-stratified permutation, B={args.B}")
        P["N2"] = label_null(eng_comb, codes, K, cands_b, obs_b, args.B, rng, patients, tools, log_every=max(args.B // 5, 1))
    if "N3" in nulls:
        log(f"N3 matched decoys, B={args.B}")
        P["N3"] = decoy_null(S_all, mat_mean, cands_b_all, obs_b, strata, decoy_pool, props_all, EXPR_PROP, args.B, rng, True, tools)
    if "N4" in nulls:
        log(f"N4 unmatched decoys, B={args.B}")
        P["N4"] = decoy_null(S_all, mat_mean, cands_b_all, obs_b, strata, decoy_pool, props_all, EXPR_PROP, args.B, rng, False, tools)
    if "N5" in nulls:
        log(f"N5 patient-replicate, B5={args.B5}")
        P["N5"] = patient_replicate_null(X, codes, K, patients, comb_idx, mat_mean, mat_max, cands_b, args.B5, rng, tools=tools)

    # ---------------- evaluations
    results = dict(meta=dict(cells=int(adata.n_obs), genes=int(adata.n_vars), types=K, patients=len(set(patients)),
                             B=args.B, B5=args.B5, tag=tag, truth_rows=int(len(truth)), truth_pos=int(truth.response.sum())),
                   validation=validation, tools={})
    mu_all = np.asarray(X.mean(0)).ravel()
    dec = pd.qcut(pd.Series(mu_all[lr_idx]).rank(method="first"), 10, labels=False).to_numpy()
    for t in tools:
        df = both[t]
        real = ~df.is_fake.to_numpy()
        n_real, n_fake = int(real.sum()), int((~real).sum())
        R = dict(n_real=n_real, n_fake=n_fake, e1={}, e2={}, e4={}, native_vs_N1=None)
        R["e1"]["ungated"] = op_metrics(op_prediction(df[real], obs_b[t][real]), truth)
        R["n_real_calls"] = {}
        if t in NATIVE_P:
            pn = liana_res[t][NATIVE_P[t]].to_numpy(dtype=float)
            R["e1"]["N6_native"] = op_metrics(op_prediction(df[real], obs_b[t][real], pn), truth)
            if "N1" in P:
                R["native_vs_N1"] = float(np.corrcoef(np.log10(pn + 1e-4), np.log10(P["N1"][t][real] + 1e-4))[0, 1])
        lig_dec = dec[pd.Index(genes[lr_idx]).get_indexer(df[real].ligand)]
        for nl in nulls:
            p = P[nl][t]
            q_real, q_fake = cn.bh(p[real]), cn.bh(p[~real])
            q_joint = cn.bh(p)
            R["e1"][nl] = op_metrics(op_prediction(df[real], obs_b[t][real], p[real]), truth)
            R.setdefault("e1q", {})[nl] = op_metrics(op_prediction(df[real], obs_b[t][real], np.where(q_real < ALPHA, 0.0, 1.0)), truth)
            R["e2"][nl] = dict(real_rate=float((q_joint[real] < ALPHA).mean()), fake_rate=float((q_joint[~real] < ALPHA).mean()),
                               real_p05=float((p[real] < 0.05).mean()), fake_p05=float((p[~real] < 0.05).mean()),
                               real_calls_sep=int((q_real < ALPHA).sum()), fake_calls_sep=int((q_fake < ALPHA).sum()))
            called = q_real < ALPHA
            R["n_real_calls"][nl] = int(called.sum())
            R["e4"][nl] = [int((called & (lig_dec == d)).sum()) for d in range(10)]
        results["tools"][t] = R
        log(f"== {t}: ungated AUPRC {R['e1']['ungated']['auprc']:.3f} | " + " | ".join(
            f"{nl} AUPRC {R['e1'][nl]['auprc']:.3f} OR {R['e1'][nl]['odds_ratio']:.2f} real {R['e2'][nl]['real_rate']:.1%} fake {R['e2'][nl]['fake_rate']:.1%}"
            for nl in nulls))
        pd.DataFrame({**{c: df[c].values for c in ["ligand_complex", "receptor_complex", "ligand", "receptor", "source", "target", "is_fake"]},
                      "score": obs_b[t], **{f"p_{nl}": P[nl][t] for nl in nulls}}).to_csv(OUT / f"rows_{t}_{tag}.csv.gz", index=False)

    # ---------------- E5 sanity: fully shuffled labels -> every null should call ~nothing
    B5s = max(args.B // 5, 50)
    if args.no_e5:
        B5s = 0
    if B5s:
        log(f"E5 shuffled labels, B={B5s}")
        S_sh = eng_comb.compute(codes_sh, K, need)
        real_c = {t: Candidates(liana_res[t], comb_pos, type_pos) for t in tools}
        real_c_all = {t: Candidates(liana_res[t], all_pos, type_pos) for t in tools}
        obs_sh = {t: tool_scores(t, S_sh, c.li, c.ri, c.si, c.ti, mat_mean) for t, c in real_c.items()}
        e5 = {"N1": label_null(eng_comb, codes_sh, K, real_c, obs_sh, B5s, rng, None, tools),
              "N2": label_null(eng_comb, codes_sh, K, real_c, obs_sh, B5s, rng, patients, tools),
              "N3": decoy_null(S_all_sh, mat_mean, real_c_all, obs_sh, strata, decoy_pool, S_all_sh["props"], EXPR_PROP, B5s, rng, True, tools),
              "N4": decoy_null(S_all_sh, mat_mean, real_c_all, obs_sh, strata, decoy_pool, S_all_sh["props"], EXPR_PROP, B5s, rng, False, tools)}
        for t in tools:
            results["tools"][t]["e5"] = {nl: dict(p05=float((e5[nl][t] < 0.05).mean()), calls=int((cn.bh(e5[nl][t]) < ALPHA).sum()))
                                         for nl in e5}
            log(f"  E5 {t}: " + " | ".join(f"{nl} p<.05 {v['p05']:.1%} calls {v['calls']}" for nl, v in results['tools'][t]['e5'].items()))

    # random baseline (Open Problems "Random Events" control)
    results["random_baseline"] = opm.random_baseline(truth, KEYS)
    (OUT / f"results_{tag}.json").write_text(json.dumps(results, indent=1, default=float))
    log("wrote", OUT / f"results_{tag}.json")


if __name__ == "__main__":
    main()
