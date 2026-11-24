"""Addendum to AMENDMENT_5 (NOT pre-registered): NATMI's Null A and Null B with B = 5,000 permutations.

Why: with B = 1,000 the smallest p is 1/1001. NATMI's BH family is ~260,000 rows (Mammary), so BH can reject
only when ~1/50 of the rows sit at that floor; Null B made no BH calls in AMENDMENT_5's power stage or on
real Mammary. Benchmark 1 met the same floor for CellChat and CellPhoneDB and ran B = 5,000 (AMENDMENT_3, A3.1).
This runs AMENDMENT_5's own code unchanged with B = 5,000. Native NATMI is not rerun: its score does not
depend on B, and its fidelity is checked by AMENDMENT_5's core.

Rows go to data/liana_nulls_b5000/ and results to benchmarks/results/liana_nulls/addendum_natmi_b5000_*.json,
so AMENDMENT_5's B = 1,000 files are never touched.

    python -m benchmarks.liana_nulls.addendum_b5000 --part core --dataset tabula_sapiens_mammary
    python -m benchmarks.liana_nulls.addendum_b5000 --part power
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from ccc_nulls import datasets

from . import amendment5, common, stages

B = 5000
NOTE = ("addendum to AMENDMENT_5, not pre-registered: AMENDMENT_5's code with B = 5,000 for Null A and Null B "
        "(p-value floor 1/5001); native NATMI not rerun")


def rows_dir() -> Path:
    d = Path(datasets.data_dir()) / "liana_nulls_b5000"
    d.mkdir(parents=True, exist_ok=True)
    return d


def save(name: str, obj: dict):
    common.save(f"addendum_natmi_b5000_{name}",
                {"amendment_5_sha256": common.amendment_hash(5), "addendum": NOTE, **obj})


def load(name: str) -> dict:
    f = common.OUT / f"addendum_natmi_b5000_{name}.json"
    return json.loads(f.read_text()) if f.exists() else {}


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--B" in argv:
        raise SystemExit("this addendum fixes B = 5,000")
    amendment5.rows_dir = rows_dir          # core rows, work dirs and downstream rows: the addendum's own folder
    stages.rows_dir = rows_dir              # biology reads core rows through stages._rows
    amendment5.save5, amendment5.load5 = save, load
    amendment5.main(argv + ["--B", str(B), "--no-native"])


if __name__ == "__main__":
    main()
