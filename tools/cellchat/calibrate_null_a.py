"""Calibration of CellChat's Null A on Tabula Sapiens with fully permuted labels (plan step 2).

    python -m tools.cellchat.calibrate_null_a --download --tissue all
    python -m tools.cellchat.calibrate_null_a --tissue mammary                  # release 2.1.2
    python -m tools.cellchat.calibrate_null_a --tissue liver --version 1.0.0     # the 2021 paper's

Shuffling every cell's label once destroys all cell-type structure, so no interaction is real.
CellChat's whole pipeline (over-expression filter, L-R selection, probabilities, Null A with its
defaults: nboot = 100, seed.use = 1) runs on each shuffled labelling, and the same with the
over-expression filter switched off. A calibrated test gives P(p < 0.05) = 5/101 and
P(p = 0) = 1/101 on the rows it reports, and BH at 5% calls nothing.

Data: two tissues of Tabula Sapiens (CELLxGENE collection e5f58829-1a66-40b5-a624-9046778e74f5,
the calibration dataset of the root README): Mammary (normal breast, closest to the TNBC atlas)
and Liver (a different organ of similar size). Raw counts -> ln(1 + CP10k), the normalisation of
the Open Problems TNBC file; cell types with < 5 cells dropped, as in LIANA.

Writes tools/cellchat/results/calibration_<tissue>_<version>.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats as sstats

from ccc_nulls import datasets
from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat import cellchat_py as cc

HERE = Path(__file__).resolve().parent
COLLECTION = datasets.TABULA_SAPIENS
# the calibration tissues, registered in ccc_nulls.datasets as tabula_sapiens_<tissue>
TISSUES = {t: {k: datasets.DATASETS[f"tabula_sapiens_{t}"][k] for k in ("title", "url", "size", "md5")}
           for t in ("mammary", "liver")}
ALPHA = 0.05


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def path(tissue: str) -> Path:
    return datasets.path(f"tabula_sapiens_{tissue}")


def download(tissue: str) -> Path:
    """``python -m ccc_nulls.datasets tabula_sapiens_<tissue>``: parallel download, md5-checked."""
    return datasets.download(f"tabula_sapiens_{tissue}")


def load(tissue: str, min_cells: int = 5):
    """``ccc_nulls.datasets.load_tabula_sapiens`` restricted to CellChatDB genes: X = ln(1 + CP10k)
    of the raw counts (library sizes use all genes), var_names = gene symbols, obs['label'] =
    cell_type. Keeping only CellChatDB genes changes nothing for CellChat, which uses
    data.signaling only, and keeps the data small."""
    return datasets.load_tabula_sapiens(f"tabula_sapiens_{tissue}",
                                        genes=cc.extract_gene(cc.subset_db(cc.load_db())), min_cells=min_cells)


def summarise(df, res, B):
    testable = df.score.to_numpy() > 0
    p = res.p_null_a_cellchat.to_numpy()[testable]          # CellChat's p: #exceed / B
    p1 = res.p_null_a.to_numpy()[testable]                  # (1 + #exceed) / (B + 1)
    k = res.n_exceed_null_a.to_numpy()[testable]
    # under H0 the observed probability's rank among the B + 1 is uniform: #exceed ~ U{0..B}
    ks = sstats.kstest((k + np.random.default_rng(0).random(len(k))) / (B + 1), "uniform") if len(k) else None
    return dict(
        n_lr=int(df.interaction_name.nunique()), n_rows=int(len(df)), n_testable=int(testable.sum()),
        frac_p_lt_0_05=float((p < ALPHA).mean()) if len(p) else None,
        frac_p_eq_0=float((p == 0).mean()) if len(p) else None,
        frac_all_rows_p_lt_0_05=float((res.p_null_a_cellchat.to_numpy() < ALPHA).mean()),
        ks_uniform_rank=dict(statistic=float(ks.statistic), pvalue=float(ks.pvalue)) if ks else None,
        calls_p_lt_0_05=int((p < ALPHA).sum()),
        calls_bh_q_lt_0_05_cellchat_p=int((bh(p) < ALPHA).sum()) if len(p) else 0,
        calls_bh_q_lt_0_05_plus_one_p=int((bh(p1) < ALPHA).sum()) if len(p) else 0,
    )


def run_once(adata, B, seed, de_filter=True, version=cc.VERSION):
    """CellChat on ``adata``. ``de_filter=False`` skips the over-expression filter: every CellChatDB
    pair whose genes are in the data is tested, so no gene is chosen with the labels it is tested on."""
    inputs = adapter.prepare_inputs(adata, "label")
    if not de_filter:
        inputs["LRsig"] = cc.identify_overexpressed_interactions(inputs["db"], inputs["genes"], inputs["genes"])
    if len(inputs["LRsig"]) == 0:
        return None
    native = adapter.run_native(inputs, nboot=2, seed=seed, version=version)   # schema rows (prob) only
    df = adapter.native_to_schema(native)
    res = nulls.null_a(adata, df, "label", B=B, seed=seed, rng="cellchat", version=version)
    return summarise(df, res, B)


def pooled(runs, B, n_boot=2000, seed=0):
    """Pooled rates over shuffles (rows weighted by shuffle), with 95% CIs from a bootstrap over
    shuffles: rows of one shuffle are correlated, shuffles are independent."""
    runs = [s for s in runs if s and s["n_testable"]]
    n = np.array([s["n_testable"] for s in runs], float)
    hit = np.array([s["frac_p_lt_0_05"] for s in runs]) * n
    zero = np.array([s["frac_p_eq_0"] for s in runs]) * n
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(runs), (n_boot, len(runs)))
    ci = lambda num: [float(x) for x in np.percentile(num[idx].sum(1) / n[idx].sum(1), [2.5, 97.5])]
    return dict(
        n_shuffles=len(runs), n_testable=int(n.sum()),
        frac_p_lt_0_05=float(hit.sum() / n.sum()), frac_p_lt_0_05_ci95=ci(hit),
        expected_p_lt_0_05=5 / (B + 1),
        frac_p_eq_0=float(zero.sum() / n.sum()), frac_p_eq_0_ci95=ci(zero), expected_p_eq_0=1 / (B + 1),
        per_shuffle_frac_p_lt_0_05=[round(s["frac_p_lt_0_05"], 4) for s in runs],
        frac_all_rows_p_lt_0_05=float(np.sum([s["frac_all_rows_p_lt_0_05"] * s["n_rows"] for s in runs])
                                      / np.sum([s["n_rows"] for s in runs])),
        calls_bh_q_lt_0_05_cellchat_p=int(sum(s["calls_bh_q_lt_0_05_cellchat_p"] for s in runs)),
        calls_bh_q_lt_0_05_plus_one_p=int(sum(s["calls_bh_q_lt_0_05_plus_one_p"] for s in runs)),
    )


def zero_hit_check(adata, version, shuffles=3, B=100, seed=1):
    """Why reported rows are anti-conservative under full permutation: split them by whether any
    permutation drove prob to exactly 0. A permuted prob of 0 always loses to an observed prob > 0,
    so conditioning on prob > 0 (what CellChat reports) favours small p. No over-expression filter."""
    from tools.cellchat.cellchat_py import r_permutations

    quantile = cc.VERSIONS[version][1]
    labels = adata.obs.label.to_numpy()
    ps, zs = [], []
    for r in range(shuffles):
        sh = adata.copy()
        sh.obs["label"] = pd.Categorical(np.random.default_rng(1000 + r).permutation(labels))
        inputs = adapter.prepare_inputs(sh, "label")
        inputs["LRsig"] = cc.identify_overexpressed_interactions(inputs["db"], inputs["genes"], inputs["genes"])
        df = adapter.native_to_schema(adapter.run_native(inputs, nboot=2, seed=seed, version=version))
        codes, levels = cc.labels_and_levels(sh.obs.label)
        K = len(levels)
        lr = nulls._lr_table(df)
        genes = pd.Index(sh.var_names)
        need = sorted({g for c in nulls.GENE_COLS for x in lr[c] for g in nulls._genes(x)} & set(genes))
        X = sparse.csc_matrix(sh.X)
        sig = genes.get_indexer(nulls.signaling_genes())
        U = sparse.csc_matrix(X[:, genes.get_indexer(need)], dtype=np.float64)
        U.data = U.data / float(X[:, sig[sig >= 0]].max())
        tm = nulls.TriMean(U, quantile=quantile)
        score = nulls.CellChatScore(lr, {g: j for j, g in enumerate(need)})
        P, perm = score(tm(codes, K)), r_permutations(len(codes), B, seed)
        exceed, zero = np.zeros(P.shape, int), np.zeros(P.shape, int)
        for b in range(B):
            Pb = score(tm(codes[perm[:, b] - 1], K))
            exceed += Pb > P
            zero += Pb == 0
        keep = P > 0
        ps.append(exceed[keep] / B)
        zs.append(zero[keep])
    p, z = np.concatenate(ps), np.concatenate(zs)
    out = {}
    for name, sel in [("all_scored_rows", np.ones_like(p, bool)), ("no_permutation_hits_0", z == 0),
                      ("some_permutations_hit_0", z > 0)]:
        out[name] = dict(n=int(sel.sum()), frac_p_lt_0_05=float(np.mean(p[sel] < ALPHA)),
                         frac_p_eq_0=float(np.mean(p[sel] == 0)))
    out["expected"] = dict(frac_p_lt_0_05=5 / (B + 1), frac_p_eq_0=1 / (B + 1))
    return out


def calibrate(tissue, version, shuffles, B, seed, zero_hit_shuffles):
    adata = load(tissue)
    out = dict(dataset=TISSUES[tissue]["title"], collection=COLLECTION, url=TISSUES[tissue]["url"],
               n_cells=int(adata.n_obs), n_genes_cellchatdb=int(adata.n_vars),
               n_groups=int(adata.obs.label.nunique()), n_donors=int(adata.obs.donor_id.nunique()),
               cellchat=cc.VERSIONS[version][0], version=version, B=B, seed=seed, shuffles=shuffles)
    log(f"{tissue}, CellChat {version}: {adata.n_obs:,} cells, {out['n_groups']} cell types, {out['n_donors']} donors")
    variants = {"cellchat_pipeline": True, "no_overexpression_filter": False}
    out["real_labels"] = {name: run_once(adata, B, seed, de, version) for name, de in variants.items()}
    log("  real labels:", {n: (s["n_testable"], round(s["frac_p_lt_0_05"], 3)) for n, s in out["real_labels"].items()})
    labels = adata.obs.label.to_numpy()
    runs = {name: [] for name in variants}
    for r in range(shuffles):
        sh = adata.copy()
        sh.obs["label"] = pd.Categorical(np.random.default_rng(1000 + r).permutation(labels))
        for name, de in variants.items():
            runs[name].append(run_once(sh, B, seed, de, version))
        log(f"  shuffle {r + 1}/{shuffles}:", {n: round(v[-1]["frac_p_lt_0_05"], 3) if v[-1] and v[-1]["n_testable"]
                                              else None for n, v in runs.items()})
    out["full_permutation"] = runs
    out["pooled"] = {name: pooled(v, B) for name, v in runs.items()}
    for name, p in out["pooled"].items():
        log(f"  pooled {name}: P(p<0.05) = {p['frac_p_lt_0_05']:.4f} {p['frac_p_lt_0_05_ci95']}, "
            f"BH calls {p['calls_bh_q_lt_0_05_cellchat_p']} / {p['calls_bh_q_lt_0_05_plus_one_p']}")
    if version != "1.0.0" and zero_hit_shuffles:
        out["zero_hit_check"] = zero_hit_check(adata, version, zero_hit_shuffles, B, seed)
        log("  zero-hit check:", {k: v.get("frac_p_lt_0_05") for k, v in out["zero_hit_check"].items()})
    f = HERE / "results" / f"calibration_{tissue}_{version}.json"
    f.write_text(json.dumps(out, indent=2))
    log("wrote", f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="mammary", choices=[*TISSUES, "all"])
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--shuffles", type=int, default=20)
    ap.add_argument("--B", type=int, default=100, help="CellChat's nboot (its default: 100)")
    ap.add_argument("--seed", type=int, default=1, help="CellChat's seed.use (its default: 1)")
    ap.add_argument("--version", default=cc.VERSION, choices=sorted(cc.VERSIONS))
    ap.add_argument("--zero-hit-shuffles", type=int, default=3)
    args = ap.parse_args()
    tissues = list(TISSUES) if args.tissue == "all" else [args.tissue]
    if args.download:
        for t in tissues:
            log("downloading", TISSUES[t]["title"])
            p = download(t)
            log("done", p, "md5", datasets.md5sum(p))
        return
    for t in tissues:
        calibrate(t, args.version, args.shuffles, args.B, args.seed, args.zero_hit_shuffles)


if __name__ == "__main__":
    main()
