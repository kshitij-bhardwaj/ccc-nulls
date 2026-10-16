"""Calibration of CellChat's Null A on Tabula Sapiens with fully permuted labels (plan step 2).

    python -m tools.cellchat.calibrate_null_a --download        # Tabula Sapiens - Mammary, 1.33 GB
    python -m tools.cellchat.calibrate_null_a                   # 5 full permutations + real labels

Shuffling every cell's label once destroys all cell-type structure, so no interaction is real.
CellChat's whole pipeline (over-expression filter, L-R selection, probabilities, Null A with its
defaults: nboot = 100, seed.use = 1) is run on each shuffled labelling. A calibrated null gives
P(p < 0.05) ~ 5/101 and P(p = 0) ~ 1/101, and BH at 5% should call ~nothing.

Data: Tabula Sapiens (CELLxGENE collection e5f58829-1a66-40b5-a624-9046778e74f5), tissue Mammary
(normal breast, the closest tissue to the TNBC evaluation atlas): raw counts -> ln(1 + CP10k), the
normalisation of the Open Problems TNBC file. Cell types with < 5 cells are dropped, as in LIANA.

Writes tools/cellchat/results/calibration_tabula_sapiens_mammary.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile

import numpy as np
from scipy import sparse
from scipy import stats as sstats

from ccc_nulls import datasets
from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat import cellchat_py as cc

HERE = Path(__file__).resolve().parent
TS_MAMMARY = dict(
    title="Tabula Sapiens - Mammary",
    url="https://datasets.cellxgene.cziscience.com/33b25f68-a9d6-4643-9283-653e937ae608.h5ad",
    size=1334205970, filename="tabula_sapiens_mammary.h5ad", groupby="cell_type",
)
ALPHA = 0.05


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def path() -> Path:
    return datasets.data_dir() / TS_MAMMARY["filename"]


def download(parts: int = 16) -> Path:
    """Parallel byte-range download (single streams from the CELLxGENE CDN can be slow)."""
    out = path()
    if out.exists() and out.stat().st_size == TS_MAMMARY["size"]:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    size = TS_MAMMARY["size"]
    step = -(-size // parts)
    with tempfile.TemporaryDirectory(dir=out.parent) as tmp:
        chunks = [(i * step, min((i + 1) * step, size) - 1, Path(tmp) / f"part{i:02d}") for i in range(parts)]
        with ThreadPoolExecutor(parts) as ex:
            list(ex.map(lambda c: datasets._fetch_range(TS_MAMMARY["url"], *c), chunks))
        with open(out, "wb") as f:
            for _, _, p in chunks:
                f.write(p.read_bytes())
    if out.stat().st_size != size:
        out.unlink()
        raise RuntimeError("size mismatch; file removed")
    return out


def _h5_strings(f, key):
    """A string column of an h5ad dataframe, plain or categorical."""
    import h5py

    o = f[key]
    if isinstance(o, h5py.Group):                          # categorical: codes + categories
        return np.asarray(o["categories"].asstr()[:], dtype=object)[o["codes"][:]]
    return np.asarray(o.asstr()[:], dtype=object)


def load(min_cells: int = 5):
    """Lean AnnData: X = ln(1 + CP10k) of the raw counts (``raw/X``), restricted to CellChatDB
    genes (library sizes use all genes), var_names = gene symbols, obs['label'] = cell_type.

    Read with h5py: the file also holds a ``scale_data`` layer too large for 8 GB of memory.
    Keeping only CellChatDB genes changes nothing for CellChat, which uses data.signaling only.
    """
    import anndata as ad
    import h5py
    import pandas as pd

    from tools.cellchat import cellchat_py as cc

    with h5py.File(path(), "r") as f:
        g = f["raw/X"]
        n_obs, n_var = len(g["indptr"]) - 1, len(f["raw/var/_index"])
        counts = sparse.csr_matrix((g["data"][:], g["indices"][:].astype(np.int32), g["indptr"][:]),
                                   shape=(n_obs, n_var))
        symbols = _h5_strings(f, "raw/var/feature_name")
        obs = pd.DataFrame({"cell_type": _h5_strings(f, "obs/cell_type"),
                            "donor_id": _h5_strings(f, "obs/donor_id")}, index=_h5_strings(f, "obs/_index"))
    lib = np.asarray(counts.sum(1, dtype=np.float64)).ravel()
    keep = np.where(pd.Index(symbols).isin(set(cc.extract_gene(cc.subset_db(cc.load_db())))))[0]
    X = sparse.csr_matrix(counts[:, keep], dtype=np.float64)
    del counts
    X = (sparse.diags(1e4 / np.where(lib > 0, lib, 1)) @ X).tocsr()
    X.data = np.log1p(X.data)
    adata = ad.AnnData(X=X.astype(np.float32), obs=obs)
    adata.var_names = pd.Index(symbols[keep]).astype(str)
    adata = adata[:, ~adata.var_names.duplicated()].copy()
    adata.obs["label"] = adata.obs["cell_type"].astype(str)
    vc = adata.obs.label.value_counts()
    adata = adata[adata.obs.label.isin(vc[vc >= min_cells].index)].copy()
    adata.obs["label"] = adata.obs.label.astype("category")
    return adata


def summarise(df, res, B):
    testable = df.score.to_numpy() > 0
    p = res.p_null_a_cellchat.to_numpy()[testable]          # CellChat's p: #exceed / B
    p1 = res.p_null_a.to_numpy()[testable]                  # (1 + #exceed) / (B + 1)
    k = res.n_exceed_null_a.to_numpy()[testable]
    # under H0 the observed probability's rank among the B + 1 is uniform: #exceed ~ U{0..B}
    ks = sstats.kstest((k + np.random.default_rng(0).random(len(k))) / (B + 1), "uniform")
    return dict(
        n_lr=int(df.interaction_name.nunique()), n_rows=int(len(df)), n_testable=int(testable.sum()),
        frac_p_lt_0_05=float((p < ALPHA).mean()), expected_p_lt_0_05=5 / (B + 1),
        frac_p_eq_0=float((p == 0).mean()), expected_p_eq_0=1 / (B + 1),
        ks_uniform_rank=dict(statistic=float(ks.statistic), pvalue=float(ks.pvalue)),
        calls_p_lt_0_05=int((p < ALPHA).sum()),
        frac_all_rows_p_lt_0_05=float((res.p_null_a_cellchat.to_numpy() < ALPHA).mean()),
        calls_bh_q_lt_0_05_cellchat_p=int((bh(p) < ALPHA).sum()),
        calls_bh_q_lt_0_05_plus_one_p=int((bh(p1) < ALPHA).sum()),
    )


def run_once(adata, B, seed, de_filter=True, version="2.2"):
    """CellChat on ``adata``. ``de_filter=False`` skips the over-expression filter: every CellChatDB
    pair whose genes are in the data is tested, so no gene is chosen with the labels it is tested on."""
    inputs = adapter.prepare_inputs(adata, "label")
    if not de_filter:
        inputs["LRsig"] = cc.identify_overexpressed_interactions(inputs["db"], inputs["genes"], inputs["genes"])
    if len(inputs["LRsig"]) == 0:
        return None
    native = adapter.run_native(inputs, nboot=2, seed=seed, version=version)    # schema rows (prob) only
    df = adapter.native_to_schema(native)
    res = nulls.null_a(adata, df, "label", B=B, seed=seed, rng="cellchat", version=version)
    return summarise(df, res, B)


def zero_hit_check(adata, shuffles=3, B=100, seed=1):
    """Why scored rows are anti-conservative under full permutation: split them by whether any
    permutation drove prob to exactly 0. A permuted prob of 0 always loses to an observed prob > 0,
    so conditioning on prob > 0 (what CellChat reports) favours small p. No over-expression filter."""
    import pandas as pd
    from scipy import sparse as sp

    from tools.cellchat.cellchat_py import r_permutations

    labels = adata.obs.label.to_numpy()
    ps, zs = [], []
    for r in range(shuffles):
        sh = adata.copy()
        sh.obs["label"] = pd.Categorical(np.random.default_rng(1000 + r).permutation(labels))
        inputs = adapter.prepare_inputs(sh, "label")
        inputs["LRsig"] = cc.identify_overexpressed_interactions(inputs["db"], inputs["genes"], inputs["genes"])
        df = adapter.native_to_schema(adapter.run_native(inputs, nboot=2, seed=seed))
        codes, levels = cc.labels_and_levels(sh.obs.label)
        K = len(levels)
        lr = nulls._lr_table(df)
        genes = pd.Index(sh.var_names)
        need = sorted({g for c in nulls.GENE_COLS for x in lr[c] for g in nulls._genes(x)} & set(genes))
        X = sp.csc_matrix(sh.X)
        sig = genes.get_indexer(nulls.signaling_genes())
        U = sp.csc_matrix(X[:, genes.get_indexer(need)], dtype=np.float64)
        U.data = U.data / float(X[:, sig[sig >= 0]].max())
        tm, score = nulls.TriMean(U), nulls.CellChatScore(lr, {g: j for j, g in enumerate(need)})
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--shuffles", type=int, default=5)
    ap.add_argument("--B", type=int, default=100, help="CellChat's nboot (its default: 100)")
    ap.add_argument("--seed", type=int, default=1, help="CellChat's seed.use (its default: 1)")
    ap.add_argument("--version", default="2.2", choices=["2.2", "1.0"], help="CellChat version of Null A")
    args = ap.parse_args()
    if args.download:
        log("downloading", TS_MAMMARY["title"])
        p = download()
        md5 = hashlib.md5()
        with open(p, "rb") as f:
            while b := f.read(1 << 24):
                md5.update(b)
        log("done", p, "md5", md5.hexdigest())
        return

    adata = load()
    out = dict(dataset=TS_MAMMARY["title"], url=TS_MAMMARY["url"], n_cells=int(adata.n_obs),
               n_genes=int(adata.n_vars), n_groups=int(adata.obs.label.nunique()),
               n_donors=int(adata.obs.donor_id.nunique()), cellchat=adapter.CELLCHAT_REF, version=args.version, B=args.B, seed=args.seed)
    log(f"{adata.n_obs:,} cells, {out['n_groups']} cell types, {out['n_donors']} donors")

    variants = {"cellchat_pipeline": True, "no_overexpression_filter": False}
    out["real_labels"] = {}
    for name, de in variants.items():
        log("real labels,", name)
        out["real_labels"][name] = s = run_once(adata, args.B, args.seed, de, args.version)
        log("  ", {k: s[k] for k in ("n_lr", "n_testable", "frac_p_lt_0_05", "calls_p_lt_0_05")})

    labels = adata.obs.label.to_numpy()
    out["full_permutation"] = {name: [] for name in variants}
    for r in range(args.shuffles):
        sh = adata.copy()
        sh.obs["label"] = np.random.default_rng(1000 + r).permutation(labels)
        sh.obs["label"] = sh.obs.label.astype("category")
        for name, de in variants.items():
            log(f"full permutation {r + 1}/{args.shuffles}, {name}")
            s = run_once(sh, args.B, args.seed, de, args.version)
            out["full_permutation"][name].append(s)
            if s:
                log("  ", {k: round(s[k], 4) if isinstance(s[k], float) else s[k] for k in (
                    "n_lr", "n_testable", "frac_p_lt_0_05", "frac_p_eq_0", "frac_all_rows_p_lt_0_05",
                    "calls_bh_q_lt_0_05_cellchat_p", "calls_bh_q_lt_0_05_plus_one_p")})

    out["pooled"] = {}
    for name in variants:
        runs = [s for s in out["full_permutation"][name] if s]
        if not runs:
            continue
        n = np.array([s["n_testable"] for s in runs])
        w = lambda key: float(np.sum([s[key] * s["n_testable"] for s in runs]) / n.sum())
        wa = lambda key: float(np.sum([s[key] * s["n_rows"] for s in runs]) / sum(s["n_rows"] for s in runs))
        out["pooled"][name] = dict(
            n_shuffles=len(runs), n_testable=int(n.sum()),
            frac_p_lt_0_05=w("frac_p_lt_0_05"), expected_p_lt_0_05=5 / (args.B + 1),
            frac_p_eq_0=w("frac_p_eq_0"), expected_p_eq_0=1 / (args.B + 1),
            frac_all_rows_p_lt_0_05=wa("frac_all_rows_p_lt_0_05"),
            per_shuffle_frac_p_lt_0_05=[round(s["frac_p_lt_0_05"], 4) for s in runs],
            calls_bh_q_lt_0_05_cellchat_p=int(sum(s["calls_bh_q_lt_0_05_cellchat_p"] for s in runs)),
            calls_bh_q_lt_0_05_plus_one_p=int(sum(s["calls_bh_q_lt_0_05_plus_one_p"] for s in runs)))
        log("pooled,", name, out["pooled"][name])
    if args.version == "2.2":
        log("zero-hit check (3 shuffles, no over-expression filter)")
        out["zero_hit_check"] = zero_hit_check(adata, B=args.B, seed=args.seed)
        log("  ", out["zero_hit_check"])
    suffix = "" if args.version == "2.2" else "_v1_0"
    f = HERE / "results" / f"calibration_tabula_sapiens_mammary{suffix}.json"
    f.write_text(json.dumps(out, indent=2))
    log("wrote", f)


if __name__ == "__main__":
    main()
