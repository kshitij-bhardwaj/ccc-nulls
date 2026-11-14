"""AMENDMENT_8.md: NATMI's specificity replaced too (variants ``swap4_<n>``), through every test of AMENDMENTS 6 and 7.

    python -m benchmarks.liana_final.amendment8 --part all
    python -m benchmarks.liana_final.amendment8 --part donor --only '^tabula_sapiens_lung\\|global\\|0$'

AMENDMENT_6's and AMENDMENT_7's parts with Connectome and NATMI captured (``lc.run(..., connectome=True,
natmi=True)``) and T2 on. Results: benchmarks/results/liana_final/amendment8/; tables: data/liana_final/amendment8/.
Every part resumes where it stopped.
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
from . import amendment4 as A4
from . import amendment6 as A6
from . import amendment7 as A7
from . import pipeline as P

OUT8 = P.OUT / "amendment8"
NULLS4, SWAP3 = A6.NULLS4, A6.SWAP3
SWAP4 = [f"swap4_{n}" for n in NULLS4]
SHOWN = [*A6.SHOWN, *SWAP4]          # liana6, swap_<n> (n != null_a), swap3_<n>, swap4_<n>


def a8_dir():
    d = P.rows_dir() / "amendment8"
    d.mkdir(parents=True, exist_ok=True)
    return d


def head() -> dict:
    return {**A7.head(), "amendment_8_sha256": hashlib.sha256((P.HERE / "AMENDMENT_8.md").read_bytes()).hexdigest()}


class Results(A2.Results):
    def __init__(self, part: str):
        self.part, self.f = part, OUT8 / f"{part}.json"
        self.obj = json.loads(self.f.read_text()) if self.f.exists() else {"runs": {}}

    def write(self, **extra):
        OUT8.mkdir(parents=True, exist_ok=True)
        self.obj.update(extra)
        h = head()
        self.f.write_text(json.dumps({**h, **{k: v for k, v in self.obj.items() if k not in h}}, indent=1, default=str))


def _run_both(*a, **kw):
    return A4._run_t2(*a, **{**kw, "connectome": True, "natmi": True})


def switch_on():
    """AMENDMENT_7's switch (AMENDMENT_6's: T2 in every run, Connectome captured, swap3 scored; A2.5's stability and
    tables), plus NATMI captured, swap4 scored, and every result and table to amendment8/."""
    A7.switch_on()
    lc.run = _run_both
    A2.Results = A6.Results = A7.Results = Results
    A6.a6_dir = a8_dir                 # AMENDMENT_6's real part and AMENDMENT_7's stability read and write here
    A6.SHOWN = A7.SHOWN = SHOWN
    A2.ALL = [v for v in A2.ALL if v not in SWAP4] + SWAP4


def table8(name: str) -> pd.DataFrame:
    return pd.read_pickle(a8_dir() / f"{name}_wide.pkl")


# ------------------------------------------------------------------ parts
def part_biology(args):
    """A2.8 for swap4 beside swap and swap3, from this amendment's real tables (AMENDMENT_6's part_biology)."""
    variants = ["swap_null_b", "swap_b2", "swap_t2", *SWAP3, *SWAP4]
    tests = {}
    for name in A2.SPATIAL:
        tbl, tr = table8(name), A2.spatial_truth(name)
        for v in variants:
            tests[(v, f"{name} spatial")] = A2.boot_or(A2.spatial_K(tbl, tr, v))
    from benchmarks.liana_nulls import downstream
    from benchmarks.liana_nulls.data import SPECIES
    for name in A2.DATASETS:
        cc = [f"call__{v}" for v in variants + [P.BASELINE]]
        calls = downstream.human_rows(table8(name)[lc.KEY + cc], SPECIES[name]) \
            .dropna(subset=["ligand_complex", "receptor_complex"]).groupby(lc.KEY, as_index=False)[cc].max()
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


def part_cases(args):
    """A2.9 for swap4: among liana6's calls (s != t), what swap4_<n> drops, and its overlap with swap3_<n>'s drops."""
    from scipy.stats import mannwhitneyu
    out = Results("cases")
    for name in A2.DATASETS:
        if name in out:
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        prop, _ = A2.type_stats(P.liana_view(adata, extra), groupby, de=False)
        del adata
        tbl = table8(name)
        m = tbl[(tbl.source != tbl.target) & tbl[f"call__{P.BASELINE}"]].reset_index(drop=True)
        B = pd.DataFrame(np.repeat(prop.ge(0.1).mean(0).to_numpy()[None], len(prop), 0), index=prop.index, columns=prop.columns)
        breadth = np.maximum(A2.per_row(m, B, "ligand_complex", "source"), A2.per_row(m, B, "receptor_complex", "target"))
        med = lambda x: float(np.median(x)) if len(x) else None
        r = {}
        for n in ("null_b", "b2", "t2"):
            drop, drop3 = ~m[f"call__swap4_{n}"].to_numpy(bool), ~m[f"call__swap3_{n}"].to_numpy(bool)
            mw = mannwhitneyu(breadth[drop], breadth[~drop], alternative="greater") if drop.any() and (~drop).any() else None
            r[f"swap4_{n}"] = dict(kept=int((~drop).sum()), dropped=int(drop.sum()), breadth_dropped_median=med(breadth[drop]),
                                   breadth_kept_median=med(breadth[~drop]), mann_whitney_p=float(mw.pvalue) if mw else None,
                                   dropped_by_both=int((drop & drop3).sum()), dropped_swap4_only=int((drop & ~drop3).sum()),
                                   dropped_swap3_only=int((~drop & drop3).sum()))
        out.add(name, r)
        log(f"cases {name}: {r}")


PARTS = {"calibration": A6.part_calibration, "lung": A6.part_lung, "power": A6.part_power,
         "hard_power": A6.part_hard_power, "mouse": A6.part_mouse, "real": A6.part_real, "biology": part_biology,
         "donor": A7.part_donor, "sensitivity": A7.part_sensitivity, "cases": part_cases,
         "tool_labels": A7.part_tool_labels, "tool_power": A7.part_tool_power, "tool_real": A7.part_tool_real}
LOCAL = ("biology", "cases")          # from the saved real tables, after every other part


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.amendment8")
    ap.add_argument("--part", choices=[*PARTS, "report", "all"], default="all")
    ap.add_argument("--only", default=None, metavar="REGEX", help="run only the runs whose tag matches (see amendment2)")
    args = ap.parse_args(argv)
    A2.SHARD = args.only
    switch_on()
    order = [p for p in PARTS if p not in LOCAL] + list(LOCAL)
    for p in (order if args.part == "all" else [args.part] if args.part != "report" else []):
        PARTS[p](args)
    if args.part in ("report", "all"):
        from .report8 import write
        write()


if __name__ == "__main__":
    main()
