"""The final pipeline: LIANA+'s consensus with our nulls in place of the tools' (plan: PLAN.md).

    python -m benchmarks.liana_final                        # everything, in order
    python -m benchmarks.liana_final --part run --dataset tnbc_wu2021
    python -m benchmarks.liana_final --part report

Parts:
  run          the five data sets: one LIANA+ run each (six methods; CellPhoneDB and CellChat captured),
               every consensus variant, F0 fidelity against LIANA+'s own rank_aggregate.
  calibration  F1: Tabula Sapiens with receptor-only or ligand-only genes shuffled, and shuffled labels.
  power        F2: planted edges on shuffled Tabula Sapiens Mammary.
  biology      F3: spatial co-location (TNBC, Allen brain, hypothalamus) and downstream response (all five);
               needs run.
  report       results/liana_final/README.md.

In code: ``ccc_nulls.liana_consensus`` (run, add_b2, variants) and ``benchmarks.liana_final.pipeline.consensus``.
Results: benchmarks/results/liana_final/*.json; tables: data/liana_final/ (not tracked).
"""
from __future__ import annotations

import argparse

from .evaluate import CAL_REPS, DATASETS, POWER_REPS, part_biology, part_calibration, part_power, part_run

PARTS = {"run": part_run, "calibration": part_calibration, "power": part_power, "biology": part_biology}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final")
    ap.add_argument("--part", choices=[*PARTS, "report", "all"], default="all")
    ap.add_argument("--dataset", action="append", dest="datasets", choices=DATASETS)
    ap.add_argument("--cal-reps", type=int, default=CAL_REPS)
    ap.add_argument("--power-reps", type=int, default=POWER_REPS)
    args = ap.parse_args(argv)
    args.datasets = args.datasets or DATASETS
    for part in (PARTS if args.part == "all" else [args.part] if args.part != "report" else []):
        PARTS[part](args)
    from .report import write
    write()


if __name__ == "__main__":
    main()
