"""Null A vs Null B in the LIANA+ pipeline, every stage, one command (plan: PLAN.md).

    python -m benchmarks.liana_nulls                       # everything, in order
    python -m benchmarks.liana_nulls --part core --dataset tnbc_wu2021
    python -m benchmarks.liana_nulls --part biology        # needs core on tnbc_wu2021 and allen_brain_atlas
    python -m benchmarks.liana_nulls --part calibration --dataset tabula_sapiens_mammary
    python -m benchmarks.liana_nulls --part power

Parts:
  core         LIANA+ (rank_aggregate and each tool's own method) + columns 1-3 per tool, written
               back into liana_res. Fidelity (T0): column 1 vs column 2, every row.
  calibration  T1 on Tabula Sapiens: known-null data (receptor or ligand genes shuffled), shuffled labels.
  power        T2: planted edges on shuffled Tabula Sapiens Mammary.
  biology      T3: spatial colocalisation (TNBC Visium, Allen brain), cytokines, fake pairs.
  amendment    Amendment 1 (AMENDMENT_1.md): hypothalamus MERFISH truth (needs core on
               moffitt2018_hypothalamus), downstream-response truths on every dataset.
  report       results/liana_nulls/README tables from the JSON files.

Needs: Python with liana 1.8, cellphonedb 5.0.1 (benchmarks/liana_nulls/requirements.txt), R with
CellChat 2.1.2 for column 1, the CellPhoneDB v5.0.0 database in data/, and the datasets of the README.
Results: benchmarks/results/liana_nulls/*.json (tracked); row tables: data/liana_nulls/ (not tracked).
"""
from __future__ import annotations

import argparse
from pathlib import Path

from ccc_nulls import datasets

from .common import AMENDMENT_DATASETS, CORE_DATASETS, cpdb_database, log, rows_dir, save


def make_tools(species: str, B: int, native: bool = True, r_cache: Path | None = None):
    from .tools import CellChat, CellPhoneDB
    return [CellChat(B=B, r_cache=r_cache, run_r=native, species=species),
            CellPhoneDB(cpdb_database(), B=B, run_native=native, species=species)]


def part_core(args):
    from . import core
    from .data import load

    for name in args.dataset or CORE_DATASETS:
        log("core:", name)
        adata, groupby, extras = load(name)
        r_cache = Path(datasets.data_dir()) / "cellchat_r" / f"{name}_nboot100_seed1"
        sp = extras["species"]
        rows, summary = core.run_dataset(adata, groupby, make_tools(sp, args.B, r_cache=r_cache), rows_dir() / name,
                                         liana_resource="consensus" if sp == "human" else f"{sp}consensus", log=log)
        for tool, df in rows.items():
            df.to_pickle(rows_dir() / f"{name}_{tool}.pkl")
        if "liana_res" in adata.uns:
            adata.uns["liana_res"].to_pickle(rows_dir() / f"{name}_liana_res.pkl")
        save(f"core_{name}", {"dataset": name, "groupby": groupby, "B": args.B, **summary},
             amendment=name in AMENDMENT_DATASETS)
        del adata, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", nargs="+", default=["core", "calibration", "power", "biology", "report"])
    ap.add_argument("--dataset", nargs="+")
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--reps", type=int, default=20, help="calibration replicates per scenario")
    ap.add_argument("--power-reps", type=int, default=5)
    args = ap.parse_args()
    parts = {"core": part_core}
    for p in args.part:
        if p in parts:
            parts[p](args)
        else:
            from . import stages
            getattr(stages, f"part_{p}")(args)


if __name__ == "__main__":
    main()
