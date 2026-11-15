"""AMENDMENT_4.md: T2 on the parts of AMENDMENT_2 that AMENDMENT_3 left out.

    python -m benchmarks.liana_final.amendment4 --part all
    python -m benchmarks.liana_final.amendment4 --part lung --only '^tabula_sapiens_lung\\|receptor\\|0$'

The runs are AMENDMENT_2's (same code, seeds, shuffles and planted pairs) with T2 switched on: ``p_t2`` is added
to every captured tool's rows right after LIANA+ runs, and every consensus that builds the Null B variants also
builds ``swap_t2`` and ``extra_t2``. Results: benchmarks/results/liana_final/amendment4/; tables:
data/liana_final/amendment4/. Every part resumes where it stopped.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time

import numpy as np
import pandas as pd

from benchmarks.liana_nulls.common import log
from ccc_nulls import liana_consensus as lc

from . import amendment2 as A2
from . import amendment3 as A3
from . import pipeline as P

OUT4 = P.OUT / "amendment4"
T2, NULL_B, SHOWN = A3.T2, A3.NULL_B, A3.SHOWN
LUNG, MAMMARY, HYPO = A2.LUNG, A2.MAMMARY, "moffitt2018_hypothalamus"


def a4_dir():
    d = P.rows_dir() / "amendment4"
    d.mkdir(parents=True, exist_ok=True)
    return d


def head() -> dict:
    return {**A3.head(), "amendment_4_sha256": hashlib.sha256((P.HERE / "AMENDMENT_4.md").read_bytes()).hexdigest()}


class Results(A2.Results):
    """AMENDMENT_2's resumable results, written to amendment4/ with AMENDMENT_4's SHA-256."""

    def __init__(self, part: str):
        self.part, self.f = part, OUT4 / f"{part}.json"
        self.obj = json.loads(self.f.read_text()) if self.f.exists() else {"runs": {}}

    def write(self, **extra):
        OUT4.mkdir(parents=True, exist_ok=True)
        self.obj.update(extra)
        h = head()
        self.f.write_text(json.dumps({**h, **{k: v for k, v in self.obj.items() if k not in h}}, indent=1, default=str))


def real_table(name: str) -> pd.DataFrame:
    """AMENDMENT_3's table of the real labels (every variant; liana6 and Null B equal the run's)."""
    return pd.read_pickle(A3.a3_dir() / f"{name}_wide.pkl")


def stability(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    """A2.stability over liana6, the Null B and the T2 variants."""
    from scipy.stats import spearmanr
    m = a.merge(b, on=lc.KEY, suffixes=("", "_2"))
    out = {"rows": int(len(m))}
    k = max(1, int(np.ceil(0.01 * len(m))))
    for v in SHOWN:
        x, y = m[f"specificity_rank__{v}"].to_numpy(), m[f"specificity_rank__{v}_2"].to_numpy()
        tx, ty = set(np.argsort(x, kind="stable")[:k]), set(np.argsort(y, kind="stable")[:k])
        out[v] = dict(spearman=float(spearmanr(x, y).statistic), top1_jaccard=len(tx & ty) / len(tx | ty))
    return out


# ------------------------------------------------------------------ T2 in every LIANA+ run
_run, _variants = lc.run, lc.variants


def _run_t2(*a, **kw):
    lrs, store = _run(*a, **kw)
    A3.add_t2(store)
    return lrs, store


def _variants_t2(lrs, store, nulls=lc.NULLS, aggregate_method="rra"):
    if "null_b" in nulls and "t2" not in nulls and aggregate_method == "rra":
        nulls = (*nulls, "t2")
    return _variants(lrs, store, nulls=nulls, aggregate_method=aggregate_method)


def switch_on():
    """T2 in every run of AMENDMENT_2's code; results to amendment4/; A2.5's stability and baseline tables
    include the T2 variants (the baseline: AMENDMENT_3's real tables)."""
    if lc.run is not _run_t2:
        lc.run, lc.variants = _run_t2, _variants_t2
    if "t2" not in lc.NULLS:
        lc.NULLS = (*lc.NULLS, "t2")
    A2.Results, A2.stability, A2.table = Results, stability, real_table
    A2.ALL = [v for v in A2.ALL if v not in T2] + T2


# ------------------------------------------------------------------ parts
def part_lung(args):
    """A2.1 on full Lung."""
    A2.calibration_part("confirm_calibration", LUNG, reps=3)
    A2.power_part("confirm_power", LUNG, effects=(1.5, 2.0), reps=5)


def part_donor(args):
    """A2.3 (a) for T2 (global permutations: T2 uses none); (b) swap_t2 and extra_t2 under within-donor permutations."""
    res = Results("donor_null")
    for name, donor in A2.DONOR.items():
        if all(f"{name}|global|{rep}" in res for rep in range(3)):
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        a0 = A2.liana_view(adata, extra)
        del adata
        for rep in range(3):
            t = f"{name}|global|{rep}"
            if t in res:
                continue
            t0 = time.time()
            a = A2.within_donor_labels(a0, groupby, donor, np.random.default_rng(9000 + rep))
            _, store = lc.run(a, groupby, resource_name=A2.resource(sp), n_perms=A2.N_PERMS, seed=A2.SEED, tools_only=True)
            r = {"dataset": name, "scheme": "global", "rep": rep, "seconds": round(time.time() - t0, 1),
                 "n_donors": int(a0.obs[donor].nunique())}
            for tool, d in store.items():
                d = d[(d.source != d.target) & (d.magnitude > 0)]
                r[tool] = {"rows": int(len(d)), "t2_untested": float(1 - d.t2_tested.mean())} | {
                    f"{n}_le_{al}": float((d[f"p_{n}"] <= al).mean()) for n in ("null_a", "null_b", "t2") for al in (0.01, 0.05)}
            res.add(t, r)
            log(f"donor_null {t}: " + ", ".join(f"{tool} A {r[tool]['null_a_le_0.05']:.3f} B {r[tool]['null_b_le_0.05']:.3f} "
                                               f"T2 {r[tool]['t2_le_0.05']:.3f}" for tool in lc.TOOLS))
        del a0
    real = Results("donor_real")
    for name in ("tnbc_wu2021", "tabula_sapiens_mammary", "tabula_sapiens_liver"):
        if name in real:
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        tbl, store, secs = A2.run_all(adata, groupby, sp, extra, strata=adata.obs[A2.DONOR[name]].astype(str),
                                      filters=False, b2=False)
        tbl.to_pickle(a4_dir() / f"donor_{name}_wide.pkl")
        m = real_table(name).merge(tbl, on=lc.KEY, suffixes=("", "_donor"))
        k = max(1, int(np.ceil(0.01 * len(m))))
        r = {"rows": int(len(m)), "seconds": secs}
        for v in ["swap_null_b", *T2]:
            tg = set(np.argsort(m[f"specificity_rank__{v}"].to_numpy(), kind="stable")[:k])
            td = set(np.argsort(m[f"specificity_rank__{v}_donor"].to_numpy(), kind="stable")[:k])
            cg, cd = m[f"call__{v}"].to_numpy(bool), m[f"call__{v}_donor"].to_numpy(bool)
            r[v] = dict(top1_jaccard=len(tg & td) / len(tg | td), calls_global=int(cg.sum()), calls_donor=int(cd.sum()),
                        call_jaccard=float((cg & cd).sum() / (cg | cd).sum()) if (cg | cd).any() else float("nan"))
        real.add(name, r)
        log(f"donor_real {name}: {r}")
        del adata, tbl, store


def part_sensitivity(args):
    """A2.5, unchanged, with T2 on."""
    A2.part_sensitivity(args)


def part_mouse(args):
    """A2.7."""
    A2.calibration_part("mouse_calibration", "allen_brain_atlas", reps=3)


def part_cases(args):
    """What swap_t2 drops among liana6's calls (s != t), its breadth, and its overlap with swap_null_b's drops."""
    from scipy.stats import mannwhitneyu
    out = Results("cases")
    for name in A2.DATASETS:
        if name in out:
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        prop, _ = A2.type_stats(A2.liana_view(adata, extra), groupby, de=False)
        del adata
        tbl = real_table(name)
        store = pd.read_pickle(P.rows_dir() / f"{name}_store.pkl")
        A3.add_t2(store)
        for tool, d in store.items():
            tbl = tbl.merge(d.drop_duplicates(lc.KEY)[lc.KEY + ["p_t2"]].rename(columns={"p_t2": f"p_t2_{tool}"}),
                            on=lc.KEY, how="left")
        m = tbl[(tbl.source != tbl.target) & tbl[f"call__{P.BASELINE}"]].reset_index(drop=True)
        B = pd.DataFrame(np.repeat(prop.ge(0.1).mean(0).to_numpy()[None], len(prop), 0), index=prop.index,
                         columns=prop.columns)
        bl, br = A2.per_row(m, B, "ligand_complex", "source"), A2.per_row(m, B, "receptor_complex", "target")
        breadth = np.maximum(bl, br)
        drop, drop_b = ~m["call__swap_t2"].to_numpy(bool), ~m["call__swap_null_b"].to_numpy(bool)
        mw = mannwhitneyu(breadth[drop], breadth[~drop], alternative="greater") if drop.any() and (~drop).any() else None
        med = lambda x: float(np.median(x)) if len(x) else None
        r = dict(kept=int((~drop).sum()), dropped=int(drop.sum()),
                 breadth_dropped_median=med(breadth[drop]), breadth_kept_median=med(breadth[~drop]),
                 breadth_ligand_dropped_median=med(bl[drop]), breadth_receptor_dropped_median=med(br[drop]),
                 mann_whitney_p=float(mw.pvalue) if mw else None,
                 dropped_by_both=int((drop & drop_b).sum()), dropped_t2_only=int((drop & ~drop_b).sum()),
                 dropped_null_b_only=int((~drop & drop_b).sum()),
                 dropped_jaccard=float((drop & drop_b).sum() / (drop | drop_b).sum()) if (drop | drop_b).any() else None)
        m = m.assign(breadth_ligand=bl.round(3), breadth_receptor=br.round(3))
        cols = lc.KEY + ["specificity_rank__liana6", "specificity_rank__swap_t2", "specificity_rank__swap_null_b",
                         *[f"p_t2_{t}" for t in lc.TOOLS], "breadth_ligand", "breadth_receptor"]
        OUT4.mkdir(parents=True, exist_ok=True)
        m[drop].sort_values("specificity_rank__liana6").head(20)[cols].to_csv(OUT4 / f"cases_{name}_dropped.csv", index=False)
        new = tbl[(tbl.source != tbl.target) & tbl["call__swap_t2"] & ~tbl[f"call__{P.BASELINE}"]]
        new.sort_values("specificity_rank__swap_t2").head(20)[[c for c in cols if not c.startswith("breadth")]] \
            .to_csv(OUT4 / f"cases_{name}_promoted.csv", index=False)
        r["promoted_rows"] = int(len(new))
        out.add(name, r)
        log(f"cases {name}: {r}")


def part_timing(args):
    """Wall time of T2 on the captured tools of the real-label run (median of 3)."""
    res = Results("timing")
    for name in (MAMMARY, HYPO):
        if name in res:
            continue
        store0 = pd.read_pickle(P.rows_dir() / f"{name}_store.pkl")
        secs = []
        for _ in range(3):
            store = {k: v.copy() for k, v in store0.items()}
            t0 = time.time()
            A3.add_t2(store)
            secs.append(time.time() - t0)
        res.add(name, dict(t2_seconds=round(float(np.median(secs)), 2), rows={t: int(len(d)) for t, d in store0.items()}))
        log(f"timing {name}: {res.runs[name]}")


PARTS = {"lung": part_lung, "donor": part_donor, "sensitivity": part_sensitivity, "mouse": part_mouse,
         "cases": part_cases, "timing": part_timing}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.amendment4")
    ap.add_argument("--part", choices=[*PARTS, "report", "all"], default="all")
    ap.add_argument("--only", default=None, metavar="REGEX", help="run only the runs whose tag matches (see amendment2)")
    args = ap.parse_args(argv)
    A2.SHARD = args.only
    switch_on()
    for p in (PARTS if args.part == "all" else [args.part] if args.part != "report" else []):
        PARTS[p](args)
    if args.part in ("report", "all"):
        from .report_amendment4 import write
        write()


if __name__ == "__main__":
    main()
