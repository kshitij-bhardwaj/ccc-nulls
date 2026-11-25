"""Rerun benchmark 1's core for CellChat or CellPhoneDB on this machine, to get their row tables back.

The row tables (data/liana_nulls/<dataset>_<tool>.pkl) are not in git; the downstream analyses need them.
This runs the same code as ``python -m benchmarks.liana_nulls --part core`` for one tool at a time, with two
differences:
- Column 1 is not rerun: CellChat 2.1.2 in R (not installed here), and CellPhoneDB's own statistical analysis
  (slow, single-threaded; CellPhoneDB runs Method 1, means only, to build its rows). Column 2, our Python
  Null A on each tool's own permutation sequence, equalled column 1 on every row of all five data sets
  (results/liana_nulls/README.md, T0).
- LIANA+'s rank_aggregate is not rerun (the row tables do not need it).
Summaries go to benchmarks/results/liana_nulls/rerun/, never over the committed core_<dataset>.json.

    CCC_NULLS_WORKERS=8 python -m benchmarks.liana_nulls.rerun_tools --tool cellchat
"""
from __future__ import annotations

import argparse
from pathlib import Path

from ccc_nulls import datasets

from . import common
from .common import AMENDMENT_DATASETS, CORE_DATASETS, cpdb_database, log, rows_dir


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_nulls.rerun_tools")
    ap.add_argument("--tool", required=True, choices=["cellchat", "cellphonedb"])
    ap.add_argument("--dataset", nargs="+", default=CORE_DATASETS + AMENDMENT_DATASETS)
    ap.add_argument("--B", type=int, default=1000)
    args = ap.parse_args(argv)
    common.OUT = common.OUT / "rerun"                     # summaries only; rows go to data/liana_nulls/ as before

    from .core import run_dataset
    from .data import decoy_universe, load
    from .tools import CellChat, CellPhoneDB

    for name in args.dataset:
        f = rows_dir() / f"{name}_{args.tool}.pkl"
        if f.exists() and (common.OUT / f"core_{name}_{args.tool}.json").exists():
            log(f"{args.tool} {name}: done")
            continue
        log(f"{args.tool}: {name}")
        adata, groupby, extras = load(name)
        decoys = decoy_universe(name)
        extra_adata = load(name, extra_genes=decoys)[0] if decoys else None
        sp = extras["species"]
        tool = (CellChat(B=args.B, run_r=False, species=sp) if args.tool == "cellchat"
                else CellPhoneDB(cpdb_database(), B=args.B, run_native=False, species=sp))
        rows, summary = run_dataset(adata, groupby, [tool], rows_dir() / name, liana_aggregate=False,
                                    liana_builtin_column=False, extra_adata=extra_adata, log=log)
        rows[args.tool].to_pickle(f)
        common.save(f"core_{name}_{args.tool}", {"dataset": name, "groupby": groupby, "B": args.B,
                                                 "rerun": "this machine; CellChat column 1 (R) not rerun", **summary},
                    amendment=name in AMENDMENT_DATASETS)
        del adata, extra_adata, rows


if __name__ == "__main__":
    main()
