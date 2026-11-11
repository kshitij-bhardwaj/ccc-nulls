"""AMENDMENT_6.md: Connectome's specificity replaced too (variants ``swap3_<n>``).

    python -m benchmarks.liana_final.amendment6 --part all
    python -m benchmarks.liana_final.amendment6 --part lung --only '^tabula_sapiens_lung\\|receptor\\|0$'

AMENDMENT_2's and AMENDMENT_3's runs with Connectome captured (``lc.run(..., connectome=True)``) and T2 on. Results:
benchmarks/results/liana_final/amendment6/; tables: data/liana_final/amendment6/. Every part resumes where it stopped.
"""
from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np
import pandas as pd

from benchmarks.liana_nulls.common import log
from ccc_nulls import liana_consensus as lc

from . import amendment2 as A2
from . import amendment3 as A3
from . import amendment4 as A4
from . import amendment5 as A5
from . import pipeline as P

OUT6 = P.OUT / "amendment6"
NULLS4 = ("null_a", "null_b", "b2", "t2")
SWAP2 = [f"swap_{n}" for n in NULLS4 if n != "null_a"]
SWAP3 = [f"swap3_{n}" for n in NULLS4]
SHOWN = ["liana6", *SWAP2, *SWAP3]


def a6_dir():
    d = P.rows_dir() / "amendment6"
    d.mkdir(parents=True, exist_ok=True)
    return d


def head() -> dict:
    return {**A5.head(), "amendment_6_sha256": hashlib.sha256((P.HERE / "AMENDMENT_6.md").read_bytes()).hexdigest()}


class Results(A2.Results):
    def __init__(self, part: str):
        self.part, self.f = part, OUT6 / f"{part}.json"
        self.obj = json.loads(self.f.read_text()) if self.f.exists() else {"runs": {}}

    def write(self, **extra):
        OUT6.mkdir(parents=True, exist_ok=True)
        self.obj.update(extra)
        h = head()
        self.f.write_text(json.dumps({**h, **{k: v for k, v in self.obj.items() if k not in h}}, indent=1, default=str))


def _run_connectome(*a, **kw):
    return A4._run_t2(*a, **{**kw, "connectome": True})


def switch_on():
    """AMENDMENT_4's switch (T2 in every run), plus Connectome captured, swap3 variants scored, results to amendment6/."""
    A4.switch_on()
    lc.run = _run_connectome
    A2.Results = Results
    A2.ALL = [v for v in A2.ALL if v not in SWAP3] + SWAP3


# ------------------------------------------------------------------ parts
def part_calibration(args):
    for name in A3.TISSUES:
        A2.calibration_part("calibration", name, reps=3, pcal=True)


def part_lung(args):
    A2.calibration_part("confirm_calibration", A2.LUNG, reps=3)
    A2.power_part("confirm_power", A2.LUNG, effects=(1.5, 2.0), reps=5)


def part_power(args):
    A2.power_part("power", A2.MAMMARY, effects=(1.5, 2.0), reps=5)


def part_hard_power(args):
    A2.power_part("hard_power", A2.MAMMARY, effects=(1.1, 1.25, 1.5), reps=5, design="hard")


def part_mouse(args):
    A2.calibration_part("mouse_calibration", "allen_brain_atlas", reps=3)


def part_real(args):
    res = Results("real")
    for name in A2.DATASETS:
        if name in res or (A2.SHARD and not __import__("re").search(A2.SHARD, name)):
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        tbl, store, secs = A2.run_all(adata, groupby, sp, extra, filters=False)
        del adata
        tbl.to_pickle(a6_dir() / f"{name}_wide.pkl")
        base = tbl[f"call__{P.BASELINE}"].to_numpy(bool)
        calls = {}
        for v in SHOWN:
            c = tbl[f"call__{v}"].to_numpy(bool)
            u = (c | base).sum()
            calls[v] = dict(calls=int(c.sum()), jaccard_with_liana6=float((c & base).sum() / u) if u else float("nan"))
        tools = {t: {f"p_{n}_lt_05": int((d[f"p_{n}"] < 0.05).sum()) for n in NULLS4} | {"rows": int(len(d))}
                 for t, d in store.items()}
        res.add(name, {"rows": int(len(tbl)), "calls": calls, "tools": tools, "seconds": secs})
        log(f"real {name}: " + ", ".join(f"{v} {calls[v]['calls']:,}" for v in SHOWN))
        del tbl, store


def part_biology(args):
    """A2.8 for the swap3 variants beside the swap variants, from the saved real tables."""
    variants = ["swap_null_b", "swap_b2", "swap_t2", *SWAP3]

    def table6(name):
        return pd.read_pickle(a6_dir() / f"{name}_wide.pkl")

    tests = {}
    for name in A2.SPATIAL:
        tbl, tr = table6(name), A2.spatial_truth(name)
        for v in variants:
            tests[(v, f"{name} spatial")] = A2.boot_or(A2.spatial_K(tbl, tr, v))
    for name in A2.DATASETS:
        calls = table6(name)[lc.KEY + [f"call__{v}" for v in variants + [P.BASELINE]]]
        for res in ("cytosig", "progeny"):
            r = pd.read_pickle(A2.rows_dir() / f"downstream_{name}_{res}.pkl")
            r = r.drop(columns=[c for c in calls if c not in lc.KEY and c in r]).merge(calls, on=lc.KEY, how="left") \
                .fillna({c: False for c in calls if c.startswith("call__")})
            for v in variants:
                tests[(v, f"{name} {res}")] = A2.boot_or(A2.downstream_K(r, v))
    ok = [k for k, t in tests.items() if t["odds_ratio"] is not None]
    for k, q in zip(ok, A2.bh(np.array([tests[k]["p"] for k in ok]))):
        tests[k]["q_all"] = float(q)
    for v in {k[0] for k in ok}:
        ks = [k for k in ok if k[0] == v]
        for k, q in zip(ks, A2.bh(np.array([tests[k]["p"] for k in ks]))):
            tests[k]["q_variant"] = float(q)
    out = Results("biology")
    out.obj["tests"] = {f"{v} | {t}": r for (v, t), r in tests.items()}
    out.write()
    log(f"biology: {len(ok)} tests; q_all < 0.05: " + ", ".join(f"{v} {t}" for (v, t) in ok if tests[(v, t)]["q_all"] < 0.05))


PARTS = {"calibration": part_calibration, "lung": part_lung, "power": part_power, "hard_power": part_hard_power,
         "mouse": part_mouse, "real": part_real, "biology": part_biology}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.amendment6")
    ap.add_argument("--part", choices=[*PARTS, "report", "all"], default="all")
    ap.add_argument("--only", default=None, metavar="REGEX", help="run only the runs whose tag matches (see amendment2)")
    args = ap.parse_args(argv)
    A2.SHARD = args.only
    switch_on()
    for p in (PARTS if args.part == "all" else [args.part] if args.part != "report" else []):
        PARTS[p](args)
    if args.part in ("report", "all"):
        from .report6 import write
        write()


if __name__ == "__main__":
    main()
