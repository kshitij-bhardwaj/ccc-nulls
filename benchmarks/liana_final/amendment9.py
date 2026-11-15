"""AMENDMENT_9.md: the false share and power at cutoffs other than the top 1%, and a calibration measure with none.

    python -m benchmarks.liana_final.amendment9 --part all
    python -m benchmarks.liana_final.amendment9 --part lung --only '^tabula_sapiens_lung\\|receptor\\|0$'

AMENDMENT_8's runs (its switch: Connectome and NATMI captured, T2 on; the same seeds, shuffles and planted pairs),
with each run scored at every cutoff in CUTS, on a curve, and by ``auroc_false``. Results:
benchmarks/results/liana_final/amendment9/, the same file names and run tags as amendment8/ (checked against them).
"""
from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np

from benchmarks.liana_nulls.common import log

from . import amendment2 as A2
from . import amendment3 as A3
from . import amendment8 as A8
from . import pipeline as P
from .evaluate import auroc

OUT9 = P.OUT / "amendment9"
CUTS = (0.005, 0.01, 0.02, 0.05, 0.10)
CURVE = np.geomspace(0.001, 1.0, 31)
_calibration_run, _power_measures = A2.calibration_run, A2.power_measures


def head() -> dict:
    return {**A8.head(), "amendment_9_sha256": hashlib.sha256((P.HERE / "AMENDMENT_9.md").read_bytes()).hexdigest()}


class Results(A2.Results):
    def __init__(self, part: str):
        self.part, self.f = part, OUT9 / f"{part}.json"
        self.obj = json.loads(self.f.read_text()) if self.f.exists() else {"runs": {}}

    def write(self, **extra):
        OUT9.mkdir(parents=True, exist_ok=True)
        self.obj.update(extra)
        h = head()
        self.f.write_text(json.dumps({**h, **{k: v for k, v in self.obj.items() if k not in h}}, indent=1, default=str))


def key(c: float) -> str:
    return f"{c:g}"


def order(tbl, v: str) -> np.ndarray:
    """Variant v's ranking, best first: evaluate.top's order (specificity_rank ascending, ties by row order)."""
    return np.argsort(tbl[f"specificity_rank__{v}"].to_numpy(), kind="stable")


def k_of(c: float, n: int) -> int:
    return max(1, int(np.ceil(c * n)))           # evaluate.top's k


def cutoff_measures(tbl, false: np.ndarray) -> dict:
    """Per variant: the false share of the top c (CUTS and CURVE) and auroc_false (no cutoff)."""
    n = len(tbl)
    out = {"false_share_at": {}, "false_share_curve": {"fractions": [float(c) for c in CURVE]}, "auroc_false": {}}
    for v in A2.variants_in(tbl):
        cum = np.cumsum(false[order(tbl, v)]) / np.arange(1, n + 1)        # false share of the top k, every k
        out["false_share_at"][v] = {key(c): float(cum[k_of(c, n) - 1]) for c in CUTS}
        out["false_share_curve"][v] = [float(cum[k_of(c, n) - 1]) for c in CURVE]
        out["auroc_false"][v] = auroc(tbl, v, ~false, false)
    return out


def calibration_run(*a, **kw) -> dict:
    """amendment2.calibration_run plus cutoff_measures, on the same table and false rows (captured from its call
    of false_rows)."""
    captured = {}

    def false_rows(tbl, scen, shuffled):
        captured["tbl"], captured["f"] = tbl, _false_rows(tbl, scen, shuffled)
        return captured["f"]

    _false_rows, A2.false_rows = A2.false_rows, false_rows
    try:
        r = _calibration_run(*a, **kw)
    finally:
        A2.false_rows = _false_rows
    return r | cutoff_measures(captured["tbl"], captured["f"])


def power_measures(tbl, plan) -> dict:
    """amendment2.power_measures plus power_at: the share of each class's planted pairs in the top c (the same
    planted rows and denominators as power_top)."""
    import pandas as pd
    out = _power_measures(tbl, plan)
    cls = plan.cls if "cls" in plan else pd.Series("single", index=plan.index)
    which = {(p.source, p.target, p.ligand, p.receptor): c for p, c in zip(plan.itertuples(), cls)}
    keys = zip(tbl.source.astype(str), tbl.target.astype(str), tbl.ligand_complex, tbl.receptor_complex)
    planted = np.array([which.get(k) for k in keys], dtype=object)
    n = len(tbl)
    for v in A2.variants_in(tbl):
        rank = np.empty(n, int)
        rank[order(tbl, v)] = np.arange(n)
        for c in sorted(set(cls)):
            m = planted == c
            out[v][c]["power_at"] = {key(q): float((rank[m] < k_of(q, n)).sum() / (cls == c).sum()) for q in CUTS}
    return out


def switch_on():
    """AMENDMENT_8's switch, then this amendment's measures and every result to amendment9/."""
    A8.switch_on()
    A2.Results = Results
    A2.calibration_run, A2.power_measures = calibration_run, power_measures


# ------------------------------------------------------------------ parts (AMENDMENT_8's, as listed in AMENDMENT_9.md)
def part_calibration(args):
    for name in A3.TISSUES:
        A2.calibration_part("calibration", name, reps=3)


def part_lung(args):
    A2.calibration_part("confirm_calibration", A2.LUNG, reps=3)


def part_mouse(args):
    A2.calibration_part("mouse_calibration", "allen_brain_atlas", reps=3)


def part_hard_power(args):
    A2.power_part("hard_power", A2.MAMMARY, effects=(1.25,), reps=5, design="hard")


PARTS = {"calibration": part_calibration, "lung": part_lung, "mouse": part_mouse, "hard_power": part_hard_power}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.amendment9")
    ap.add_argument("--part", choices=[*PARTS, "report", "all"], default="all")
    ap.add_argument("--only", default=None, metavar="REGEX", help="run only the runs whose tag matches (see amendment2)")
    args = ap.parse_args(argv)
    A2.SHARD = args.only
    switch_on()
    for p in (list(PARTS) if args.part == "all" else [args.part] if args.part != "report" else []):
        PARTS[p](args)
        log(f"amendment9 {p}: done")
    if args.part in ("report", "all"):
        from .report9 import write
        write()


if __name__ == "__main__":
    main()
