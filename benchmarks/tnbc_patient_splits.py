"""E3: are the significant calls reproducible across independent patients?

Split the 10 TNBC patients into two disjoint halves (5 vs 5), run every null on each half
(same candidate rows, the ones LIANA kept on the full data), and compare the halves:
  - Spearman correlation of -log10 p across interactions observed in both halves;
  - Jaccard overlap of each half's top 5% interactions (ranked by p, ties by tool score);
  - Jaccard of BH calls at q < 0.05 (undefined when a null makes no calls).
Repeated for several random splits.

A null whose calls reflect real, shared biology should give overlapping sets; calls driven
by noise or by one patient's composition should not replicate.

python benchmarks/tnbc_patient_splits.py --splits 3 --B 500
"""
import argparse
import json
import pickle
import time
import warnings
from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import ccc_nulls as cn
from ccc_nulls.multitool import TOOLS, Candidates, StatsEngine, decoy_null, label_null, tool_scores

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "tnbc_wu2021_log_cp10k.h5ad"
OUT = ROOT / "benchmarks" / "results" / "tnbc"
EXPR_PROP, ALPHA = 0.10, 0.05


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def jaccard(a, b):
    u = (a | b).sum()
    return float((a & b).sum() / u) if u else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", type=int, default=3)
    ap.add_argument("--B", type=int, default=500)
    ap.add_argument("--tools", default="CellPhoneDB,Connectome,NATMI,SingleCellSignalR,log2FC")
    args = ap.parse_args()
    tools = args.tools.split(",")
    rng = np.random.default_rng(11)

    raw = ad.read_h5ad(DATA)
    X_full = raw.layers["normalized"].astype(np.float32).tocsr()
    obs = raw.obs.copy()
    genes = pd.Index(raw.var_names.astype(str))
    del raw
    obs["label"] = obs["label"].astype(str)
    patients = obs["orig.ident"].astype(str).to_numpy()
    pats = np.unique(patients)
    liana_res = pickle.loads((OUT / "liana_full.pkl").read_bytes())   # candidate rows from the full run
    resource = li.resource.select_resource("consensus")
    res_genes = set(resource.ligand.str.split("_").explode()) | set(resource.receptor.str.split("_").explode())
    types = sorted(obs.label.unique())
    type_pos = {t: i for i, t in enumerate(types)}
    all_pos = {g: i for i, g in enumerate(genes)}

    out = {t: {nl: [] for nl in ["N1", "N2", "N3", "N4"]} for t in tools}
    for sp in range(args.splits):
        perm = rng.permutation(pats)
        halves = [perm[:5], perm[5:]]
        calls = {t: {nl: [] for nl in out[t]} for t in tools}
        for h, hp in enumerate(halves):
            cells = np.isin(patients, hp)
            X = X_full[cells]
            lab = obs.label.to_numpy()[cells]
            codes = np.array([type_pos[x] for x in lab])
            K = len(types)
            present = np.bincount(codes, minlength=K) >= 5
            nonempty = int((np.asarray(X.sum(0)).ravel() > 0).sum())
            mat_mean = float(X.sum(dtype="float64") / (X.shape[0] * nonempty))
            mat_max = float(X.max())
            lr_idx = np.where(genes.isin(list(res_genes)))[0]
            lr_pos = {g: i for i, g in enumerate(genes[lr_idx])}
            need = set().union(*(TOOLS[t][2] for t in tools)) | {"props"}
            eng = StatsEngine(X, lr_idx, mat_mean, mat_max)
            eng_all = StatsEngine(X, np.arange(len(genes)), mat_mean, mat_max)
            S = eng.compute(codes, K, need)
            S_all = eng_all.compute(codes, K, need)
            strata = cn.gene_strata(X)
            pool = np.ones(len(genes), bool)
            pool[[all_pos[g] for g in res_genes if g in all_pos]] = False
            cands = {t: Candidates(liana_res[t], lr_pos, type_pos) for t in tools}
            cands_all = {t: Candidates(liana_res[t], all_pos, type_pos) for t in tools}
            obs_s = {t: tool_scores(t, S, c.li, c.ri, c.si, c.ti, mat_mean) for t, c in cands.items()}
            log(f"split {sp} half {h}: {X.shape[0]:,} cells, patients {list(hp)}")
            P = {"N1": label_null(eng, codes, K, cands, obs_s, args.B, rng, None, tools),
                 "N2": label_null(eng, codes, K, cands, obs_s, args.B, rng, patients[cells], tools),
                 "N3": decoy_null(S_all, mat_mean, cands_all, obs_s, strata, pool, S_all["props"], EXPR_PROP, args.B, rng, True, tools),
                 "N4": decoy_null(S_all, mat_mean, cands_all, obs_s, strata, pool, S_all["props"], EXPR_PROP, args.B, rng, False, tools)}
            for t in tools:
                c = cands[t]
                ok = present[c.si] & present[c.ti]           # both types observed in this half
                for nl in P:
                    p = np.where(ok, P[nl][t], 1.0)
                    calls[t][nl].append((p, obs_s[t], ok))
        for t in tools:
            for nl in calls[t]:
                (pa, sa, oka), (pb, sb, okb) = calls[t][nl]
                both = oka & okb
                rho = float(spearmanr(-np.log10(pa[both]), -np.log10(pb[both])).correlation)
                k = int(0.05 * both.sum())

                def top(p, s):
                    order = np.lexsort((-s[both], p[both]))       # p ascending, then score descending
                    m = np.zeros(both.sum(), bool); m[order[:k]] = True
                    return m
                ta, tb = top(pa, sa), top(pb, sb)
                ca, cb = cn.bh(pa) < ALPHA, cn.bh(pb) < ALPHA
                out[t][nl].append(dict(spearman=rho, top5_jaccard=jaccard(ta, tb), n_both=int(both.sum()),
                                       bh_jaccard=jaccard(ca, cb), n_a=int(ca.sum()), n_b=int(cb.sum())))
            log(f"  split {sp} {t}: " + " | ".join(
                f"{nl} rho={out[t][nl][-1]['spearman']:.2f} top5J={out[t][nl][-1]['top5_jaccard']:.2f}" for nl in out[t]))
    clean = {t: {nl: [{k: (None if isinstance(v, float) and v != v else v) for k, v in r.items()} for r in rs]
                 for nl, rs in d.items()} for t, d in out.items()}
    (OUT / "patient_splits.json").write_text(json.dumps(dict(B=args.B, splits=args.splits, results=clean), indent=1))
    log("wrote", OUT / "patient_splits.json")


if __name__ == "__main__":
    main()
