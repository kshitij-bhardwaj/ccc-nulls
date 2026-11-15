"""AMENDMENT_3.md: T2 (PR #11), the restandardised cell-type-ID null, inside LIANA+'s consensus.

    python -m benchmarks.liana_final.amendment3 --part all
    python -m benchmarks.liana_final.amendment3 --part calibration --only '^tabula_sapiens_mammary\\|'

T2 replaces the p-value only, like PLAN.md's nulls: ``p_t2`` is added to each captured tool's rows, and
``swap_t2`` / ``extra_t2`` are built by ``liana_consensus.variants`` like ``swap_null_b`` / ``extra_null_b``.
The runs are AMENDMENT_2's (same seeds, shuffles and planted pairs) with T2 switched on, so liana6 and the
Null B variants come out of the same runs. Results: benchmarks/results/liana_final/amendment3/; tables:
data/liana_final/amendment3/. Every part resumes where it stopped.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time

import numpy as np
import pandas as pd

from benchmarks.liana_nulls.common import log
from ccc_nulls import liana_consensus as lc

from . import amendment2 as A2
from . import evaluate as E
from . import pipeline as P

OUT3 = P.OUT / "amendment3"
T2 = ["swap_t2", "extra_t2"]
NULL_B = ["swap_null_b", "extra_null_b"]
SHOWN = ["liana6", *NULL_B, *T2]
MIN_REF = 3
TISSUES = ["tabula_sapiens_mammary", "tabula_sapiens_liver"]


def a3_dir():
    d = P.rows_dir() / "amendment3"
    d.mkdir(parents=True, exist_ok=True)
    return d


def head() -> dict:
    sha = lambda n: hashlib.sha256((P.HERE / f"AMENDMENT_{n}.md").read_bytes()).hexdigest()
    return {"plan_sha256": P.plan_hash(), "amendment_1_sha256": sha(1), "amendment_2_sha256": sha(2),
            "amendment_3_sha256": sha(3), "written": time.strftime("%Y-%m-%d %H:%M")}


# ------------------------------------------------------------------ T2 on the captured tools
def add_t2(store: dict, min_ref: int = MIN_REF) -> None:
    """``p_t2`` for every captured row: T2 on the tool's magnitude over the interaction's other scored
    (source, target) rows; rows with fewer than ``min_ref`` references get p = 1 (``t2_tested`` False)."""
    from ccc_nulls.nullb_type_restd import p_values
    for tool, d in store.items():
        u = d.drop_duplicates(lc.KEY).reset_index(drop=True)
        p = p_values(u.assign(tested=True), score="magnitude", min_ref=min_ref)
        t2 = u[lc.KEY].assign(p_t2=np.where(np.isfinite(p), p, 1.0), t2_tested=np.isfinite(p))
        store[tool] = d.drop(columns=[c for c in ("p_t2", "t2_tested") if c in d]).merge(t2, on=lc.KEY, how="left")


_add_b2 = lc.add_b2


def _add_b2_and_t2(store, adata, groupby, excluded, seed=0):
    _add_b2(store, adata, groupby, excluded, seed)
    add_t2(store)


class Results(A2.Results):
    """AMENDMENT_2's resumable results, written to amendment3/ with AMENDMENT_3's SHA-256."""

    def __init__(self, part: str):
        self.part, self.f = part, OUT3 / f"{part}.json"
        self.obj = json.loads(self.f.read_text()) if self.f.exists() else {"runs": {}}

    def write(self, **extra):
        OUT3.mkdir(parents=True, exist_ok=True)
        self.obj.update(extra)
        h = head()
        self.f.write_text(json.dumps({**h, **{k: v for k, v in self.obj.items() if k not in h}}, indent=1, default=str))


def save3(name: str, obj: dict):
    OUT3.mkdir(parents=True, exist_ok=True)
    (OUT3 / f"{name}.json").write_text(json.dumps({**head(), **obj}, indent=1, default=str))
    log("wrote", OUT3 / f"{name}.json")


def switch_on():
    """T2 in every run of AMENDMENT_2's code: lc.NULLS gains "t2" (so run_all builds swap_t2 and extra_t2 and
    p_calibration covers p_t2), add_b2 also adds p_t2, results go to amendment3/."""
    if "t2" not in lc.NULLS:
        lc.NULLS = (*lc.NULLS, "t2")
        lc.add_b2 = _add_b2_and_t2
    A2.Results = Results
    A2.ALL = [v for v in A2.ALL if v not in T2] + T2


# ------------------------------------------------------------------ parts
def part_calibration(args):
    """F1 + A2.6: Mammary and Liver, both sides, 3 replicates, tool-level p kept."""
    for name in TISSUES:
        A2.calibration_part("calibration", name, reps=3, pcal=True)


def part_power(args):
    """F2's design (power.background / design / plant, 10,000 cells), f = 1.5 and 2.0, 5 replicates."""
    A2.power_part("power", A2.MAMMARY, effects=(1.5, 2.0), reps=5)


def part_hard_power(args):
    """A2.2's design."""
    A2.power_part("hard_power", A2.MAMMARY, effects=(1.1, 1.25, 1.5), reps=5, design="hard")


def part_real(args):
    """The five data sets: every variant, saved for biology; the earlier variants checked against the run tables."""
    res = Results("real")
    for name in E.DATASETS:
        if name in res:
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        tbl, store, secs = A2.run_all(adata, groupby, sp, extra, filters=False)
        del adata
        tbl.to_pickle(a3_dir() / f"{name}_wide.pkl")
        base = tbl[f"call__{P.BASELINE}"].to_numpy(bool)
        calls = {}
        for v in SHOWN:
            c = tbl[f"call__{v}"].to_numpy(bool)
            u = (c | base).sum()
            calls[v] = dict(calls=int(c.sum()), jaccard_with_liana6=float((c & base).sum() / u) if u else float("nan"))
        tools = {t: dict(rows=int(len(d)), t2_untested=float(1 - d.t2_tested.mean()), p_t2_lt_05=int((d.p_t2 < 0.05).sum()),
                         p_null_b_lt_05=int((d.p_null_b < 0.05).sum()), p_null_a_lt_05=int((d.p_null_a < 0.05).sum()))
                 for t, d in store.items()}
        rec = {"rows": int(len(tbl)), "calls": calls, "tools": tools, "seconds": secs}
        old_f = P.rows_dir() / f"{name}_wide.pkl"
        if old_f.exists():
            old = pd.read_pickle(old_f)
            m = old.merge(tbl, on=lc.KEY, suffixes=("_old", "_new"))
            rec["same_as_run"] = {v: bool(len(m) == len(old) == len(tbl) and np.array_equal(
                m[f"specificity_rank__{v}_old"].to_numpy(), m[f"specificity_rank__{v}_new"].to_numpy())) for v in ["liana6", *NULL_B]}
        res.add(name, rec)
        log(f"real {name}: " + ", ".join(f"{v} {calls[v]['calls']:,}" for v in SHOWN) + f"; same as run: {rec.get('same_as_run')}")
        del tbl, store


def part_biology(args):
    """F3 (evaluate.part_biology unchanged) on the tables of ``real``, with the T2 variants added."""
    E.rows_dir, E.VARIANTS, E.save = a3_dir, [*P.VARIANTS, *T2], save3
    E.part_biology(args)


# ------------------------------------------------------------------ decision and head-to-head
def _load(n):
    f = OUT3 / f"{n}.json"
    return json.loads(f.read_text()) if f.exists() else None


def decision() -> dict:
    """AMENDMENT_1's rule for swap_t2 and extra_t2 (None where a part is missing)."""
    from benchmarks.liana_nulls.stages import mean_ci
    cal, pw, bio = _load("calibration"), _load("power"), _load("biology")
    out = {}
    for v in T2:
        c1 = None
        if cal:
            cells = [(t, s) for t in TISSUES for s in ("receptor", "ligand")]
            runs = lambda t, s: [r for r in cal["runs"].values() if r["dataset"] == t and r["scenario"] == s]
            if all(len(runs(t, s)) == 3 for t, s in cells):
                c1 = all(np.mean([r["false_share_top"][v] for r in runs(t, s)]) <
                         np.mean([r["false_share_top"][P.BASELINE] for r in runs(t, s)]) for t, s in cells)
        c2 = None
        if pw:
            rs = [r for r in pw["runs"].values() if r["effect"] == 2.0]
            if len(rs) == 5:
                c2 = mean_ci([r[v]["single"]["power_top"] for r in rs])[0] >= 0.8 and \
                    mean_ci([r[v]["single"]["auroc"] for r in rs])[0] >= mean_ci([r[P.BASELINE]["single"]["auroc"] for r in rs])[0]
        c3 = None
        if bio:
            cis = [bio["spatial"][d][v]["ci95"] for d in bio["spatial"] if bio["spatial"][d][v].get("ci95")]
            c3 = any(lo > 1 for lo, hi in cis) and not any(hi < 1 for lo, hi in cis)
        verdict = None if None in (c1, c2, c3) else ("improves LIANA+" if c1 and c2 and c3 else
                                                      "valid, not biologically distinguishable" if c1 and c2 else "not supported")
        out[v] = dict(calibration=c1, power=c2, biology=c3, verdict=verdict)
    return out


def part_report(args):
    from .report_amendment3 import write
    write()


PARTS = {"calibration": part_calibration, "power": part_power, "hard_power": part_hard_power, "real": part_real,
         "biology": part_biology}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.amendment3")
    ap.add_argument("--part", choices=[*PARTS, "report", "all"], default="all")
    ap.add_argument("--only", default=None, metavar="REGEX", help="run only the runs whose tag matches (see amendment2)")
    args = ap.parse_args(argv)
    A2.SHARD = args.only
    switch_on()
    for p in (PARTS if args.part == "all" else [args.part] if args.part != "report" else []):
        PARTS[p](args)
    if args.part in ("report", "all"):
        part_report(args)


if __name__ == "__main__":
    main()
