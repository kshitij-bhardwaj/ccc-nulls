"""Run one tool's Null A / Null B on an Open Problems dataset and score them the same way for every tool.

    python benchmarks/run_tool.py --tool cellphonedb                       # LIANA version, TNBC
    python benchmarks/run_tool.py --tool cellchat --source native --B 100
    python benchmarks/run_tool.py --tool natmi --quick --B 50              # 8,000-cell smoke test

Steps: load the dataset -> interactions from tools/<tool>/adapter.py (LIANA reference or native)
-> tools/<tool>/nulls.py null_a / null_b (whichever are implemented) -> ccc_nulls.schema checks
-> summary. Writes to benchmarks/results/tools/:
    <tool>_<source>_<tag>.json       summary (calls, agreement with the native p, Open Problems metrics)
    <tool>_<source>_<tag>.csv.gz     one row per interaction with every p-value (git-ignored)
"""
import argparse
import importlib
import json
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from scipy import stats as sstats

from ccc_nulls import datasets, schema
from ccc_nulls.evaluation import openproblems as opm
from ccc_nulls.stats import bh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))                    # so tools/<tool> imports when run as a script
OUT =ROOT / "benchmarks" / "results" / "tools"
TOOLS = ["cellphonedb", "cellchat", "natmi", "connectome"]
ALPHA = 0.05


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def finite(x):
    return None if x is None or not np.isfinite(x) else float(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True, choices=TOOLS)
    ap.add_argument("--source", default="liana", choices=["liana", "native"])
    ap.add_argument("--dataset", default="tnbc_wu2021", choices=["tnbc_wu2021", "allen_brain_atlas"])
    ap.add_argument("--B", type=int, default=1000, help="permutations / decoy draws for each null")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--quick", action="store_true", help="subsample 8,000 cells")
    ap.add_argument("--nulls", default="a,b")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()
    tag = args.tag or ("quick" if args.quick else "full")
    adapter = importlib.import_module(f"tools.{args.tool}.adapter")
    nulls = importlib.import_module(f"tools.{args.tool}.nulls")
    groupby = datasets.DATASETS[args.dataset]["groupby"]
    keys = datasets.DATASETS[args.dataset]["truth_keys"]

    log("loading", args.dataset)
    adata, truth = datasets.load_openproblems(args.dataset)
    if args.quick:
        rng = np.random.default_rng(args.seed)
        adata = adata[np.sort(rng.choice(adata.n_obs, min(8000, adata.n_obs), replace=False))].copy()
        vc = adata.obs[groupby].value_counts()
        adata = adata[adata.obs[groupby].isin(vc[vc >= 5].index)].copy()
        adata.obs[groupby] = adata.obs[groupby].astype(str).astype("category")
    log(f"{adata.n_obs:,} cells x {adata.n_vars:,} genes, {adata.obs[groupby].nunique()} cell types")

    log("interactions:", args.tool, args.source)
    if args.source == "liana":
        df = adapter.liana_reference(adata, groupby, seed=args.seed)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            df = adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(adata, groupby, Path(tmp)),
                                                             Path(tmp)))
    df = schema.validate(df)
    log(f"  {len(df):,} interactions")

    pcols = {}
    for which in args.nulls.split(","):
        col = schema.NULL_P[which]
        t0 = time.time()
        try:
            res = getattr(nulls, f"null_{which}")(adata, df.copy(), groupby, B=args.B, seed=args.seed)
        except NotImplementedError as e:
            log(f"  null {which.upper()}: skipped ({e})")
            continue
        df[col] = schema.validate_null(df, res, which)[col].to_numpy()
        pcols[f"null_{which}"] = col
        log(f"  null {which.upper()}: {time.time() - t0:.0f} s")
    if df["p_native"].notna().any():
        pcols["native"] = "p_native"

    summary = dict(tool=args.tool, source=args.source, dataset=args.dataset, B=args.B, seed=args.seed,
                   n_cells=int(adata.n_obs), n_interactions=int(len(df)), score_name=str(df.score_name.iloc[0]),
                   calls={}, vs_native={}, openproblems={})
    for name, col in pcols.items():
        p = df[col].to_numpy(dtype=float)
        summary["calls"][name] = dict(p_lt_0_05=int((p < ALPHA).sum()), bh_q_lt_0_05=int((bh(p) < ALPHA).sum()))
    if "native" in pcols:
        pn = df.p_native.to_numpy(dtype=float)
        for name, col in pcols.items():
            if name == "native":
                continue
            p = df[col].to_numpy(dtype=float)
            a, b = p < ALPHA, pn < ALPHA
            summary["vs_native"][name] = dict(
                spearman=finite(sstats.spearmanr(p, pn).statistic), both=int((a & b).sum()),
                only_ours=int((a & ~b).sum()), only_native=int((~a & b).sum()))

    # Open Problems protocol: the tool score, kept only where p <= 0.05 when a p-value is used
    score = df.score.to_numpy(dtype=float)
    summary["openproblems"]["score_only"] = opm.metrics(opm.predictions(df, score, keys), truth, keys)
    for name, col in pcols.items():
        summary["openproblems"][name] = opm.metrics(opm.predictions(df, score, keys, p=df[col].to_numpy()), truth, keys)
    summary["openproblems"]["random"] = opm.random_baseline(truth, keys)
    summary["openproblems"] = {k: {m: finite(v) for m, v in d.items()} for k, d in summary["openproblems"].items()}

    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / f"{args.tool}_{args.source}_{tag}"
    df.to_csv(f"{stem}.csv.gz", index=False)
    Path(f"{stem}.json").write_text(json.dumps(summary, indent=2))
    log("wrote", f"{stem}.json")
    print(json.dumps({k: summary[k] for k in ("calls", "vs_native", "openproblems")}, indent=2))


if __name__ == "__main__":
    main()
