"""AMENDMENT_5.md: the tests Null B had and T2 or B2 missed.

    python -m benchmarks.liana_final.amendment5 [--part biology|donor_b2|cases_b2|all]

Results: benchmarks/results/liana_final/amendment5/, each file with AMENDMENT_5's SHA-256. Every part resumes where
it stopped.
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
from . import amendment4 as A4
from . import pipeline as P

OUT5 = P.OUT / "amendment5"
T2 = A3.T2
N_PERMS_SCORE = 10          # B2 uses no permutations; LIANA+ needs some to run, and they change no score


def head() -> dict:
    return {**A4.head(), "amendment_5_sha256": hashlib.sha256((P.HERE / "AMENDMENT_5.md").read_bytes()).hexdigest()}


class Results(A2.Results):
    def __init__(self, part: str):
        self.part, self.f = part, OUT5 / f"{part}.json"
        self.obj = json.loads(self.f.read_text()) if self.f.exists() else {"runs": {}}

    def write(self, **extra):
        OUT5.mkdir(parents=True, exist_ok=True)
        self.obj.update(extra)
        h = head()
        self.f.write_text(json.dumps({**h, **{k: v for k, v in self.obj.items() if k not in h}}, indent=1, default=str))


def with_t2(name: str, tbl: pd.DataFrame) -> pd.DataFrame:
    """``tbl`` with AMENDMENT_3's T2 calls for the same rows."""
    t2 = A4.real_table(name)[lc.KEY + [f"call__{v}" for v in T2]]
    return tbl.drop(columns=[c for c in t2 if c not in lc.KEY and c in tbl]).merge(t2, on=lc.KEY, how="left") \
        .fillna({f"call__{v}": False for v in T2})


def donor_calls() -> dict:
    """TNBC's donor-aware variants: swap_null_b (AMENDMENT_2) and extra_t2 (AMENDMENT_4)."""
    out = {}
    for d, col, new in ((A2.a2_dir(), "swap_null_b", "swap_null_b_donor"), (A4.a4_dir(), "extra_t2", "extra_t2_donor")):
        f = d / "donor_tnbc_wu2021_wide.pkl"
        if f.exists():
            out[new] = pd.read_pickle(f)[lc.KEY + [f"call__{col}"]].rename(columns={f"call__{col}": f"call__{new}"})
    return out


# ------------------------------------------------------------------ A2.8 with T2
def part_biology(args):
    variants = [v for v in A2.VARIANTS if v != A2.BASELINE] + T2
    donors = donor_calls()

    def add_donor(t, name):
        if name != "tnbc_wu2021":
            return t, []
        for new, d in donors.items():
            t = t.merge(d, on=lc.KEY, how="left").fillna({f"call__{new}": False})
        return t, list(donors)

    tests = {}
    for name in A2.SPATIAL:
        tbl, extra = add_donor(with_t2(name, A2.table(name)), name)
        tr = A2.spatial_truth(name)
        for v in variants + extra:
            tests[(v, f"{name} spatial")] = A2.boot_or(A2.spatial_K(tbl, tr, v))
    for name in A2.DATASETS:
        for res in ("cytosig", "progeny"):
            r, extra = add_donor(with_t2(name, pd.read_pickle(A2.rows_dir() / f"downstream_{name}_{res}.pkl")), name)
            for v in variants + extra:
                tests[(v, f"{name} {res}")] = A2.boot_or(A2.downstream_K(r, v))
    ok = [k for k, t in tests.items() if t["odds_ratio"] is not None]
    for k, q in zip(ok, A2.bh(np.array([tests[k]["p"] for k in ok]))):
        tests[k]["q_all"] = float(q)
    for v in {k[0] for k in ok}:
        ks = [k for k in ok if k[0] == v]
        for k, q in zip(ks, A2.bh(np.array([tests[k]["p"] for k in ks]))):
            tests[k]["q_variant"] = float(q)
    # check: AMENDMENT_2's odds ratios and p-values are reproduced
    a2 = json.loads((A2.OUT2 / "biology.json").read_text())["tests"]
    same = [k for k in a2 if k in {f"{v} | {t}" for v, t in ok}
            and all(a2[k].get(f) == tests[tuple(k.split(" | "))].get(f) for f in ("odds_ratio", "p"))]
    out = Results("biology")
    out.obj["tests"] = {f"{v} | {t}": r for (v, t), r in tests.items()}
    out.write(check_amendment2=dict(identical=len(same), compared=len([k for k in a2 if a2[k].get("odds_ratio") is not None])))
    log(f"biology: {len(ok)} tests; AMENDMENT_2 reproduced in {len(same)}; q_all < 0.05: "
        + ", ".join(f"{v} {t}" for (v, t) in ok if tests[(v, t)]["q_all"] < 0.05))


# ------------------------------------------------------------------ A2.3 (a) for B2
def part_donor_b2(args):
    from benchmarks.liana_nulls.extra import excluded_genes
    res = Results("donor_b2")
    for name, donor in A2.DONOR.items():
        if getattr(args, "dataset", None) and name != args.dataset:
            continue
        if all(f"{name}|{rep}" in res for rep in range(3)):
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        a0 = A2.liana_view(adata, extra)
        for rep in range(3):
            t = f"{name}|{rep}"
            if t in res:
                continue
            t0 = time.time()
            a = A2.within_donor_labels(a0, groupby, donor, np.random.default_rng(9000 + rep))
            full = adata.copy()
            full.obs[groupby] = a.obs[groupby].to_numpy()
            _, store = lc.run(a, groupby, resource_name=A2.resource(sp), n_perms=N_PERMS_SCORE, seed=A2.SEED, tools_only=True)
            lc.add_b2(store, full, groupby, excluded_genes())
            del full
            r = {"dataset": name, "rep": rep, "seconds": round(time.time() - t0, 1), "n_donors": int(a0.obs[donor].nunique())}
            for tool, d in store.items():
                d = d[(d.source != d.target) & (d.magnitude > 0)]
                p = d["p_b2"].to_numpy(float)
                tested = np.isfinite(p)
                r[tool] = {"rows": int(len(d)), "rows_tested": int(tested.sum())} | {
                    f"b2_le_{al}": float((np.nan_to_num(p, nan=1.0) <= al).mean()) for al in (0.01, 0.05)} | {
                    f"b2_tested_le_{al}": float((p[tested] <= al).mean()) if tested.any() else None for al in (0.01, 0.05)}
            res.add(t, r)
            log(f"donor_b2 {t}: " + ", ".join(f"{tool} {r[tool]['b2_le_0.05']:.3f} (tested {r[tool]['b2_tested_le_0.05']})"
                                              for tool in lc.TOOLS))
        del adata, a0


# ------------------------------------------------------------------ A2.9 for B2
def part_cases_b2(args):
    from scipy.stats import mannwhitneyu
    out = Results("cases_b2")
    for name in A2.DATASETS:
        if name in out:
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        prop, _ = A2.type_stats(A2.liana_view(adata, extra), groupby, de=False)
        del adata
        tbl = A4.real_table(name)
        m = tbl[(tbl.source != tbl.target) & tbl[f"call__{P.BASELINE}"]].reset_index(drop=True)
        B = pd.DataFrame(np.repeat(prop.ge(0.1).mean(0).to_numpy()[None], len(prop), 0), index=prop.index, columns=prop.columns)
        breadth = np.maximum(A2.per_row(m, B, "ligand_complex", "source"), A2.per_row(m, B, "receptor_complex", "target"))
        drop, drop_b = ~m["call__swap_b2"].to_numpy(bool), ~m["call__swap_null_b"].to_numpy(bool)
        mw = mannwhitneyu(breadth[drop], breadth[~drop], alternative="greater") if drop.any() and (~drop).any() else None
        med = lambda x: float(np.median(x)) if len(x) else None
        out.add(name, dict(kept=int((~drop).sum()), dropped=int(drop.sum()),
                           breadth_dropped_median=med(breadth[drop]), breadth_kept_median=med(breadth[~drop]),
                           mann_whitney_p=float(mw.pvalue) if mw else None,
                           dropped_by_both=int((drop & drop_b).sum()), dropped_b2_only=int((drop & ~drop_b).sum()),
                           dropped_null_b_only=int((~drop & drop_b).sum()),
                           dropped_jaccard=float((drop & drop_b).sum() / (drop | drop_b).sum()) if (drop | drop_b).any() else None))
        log(f"cases_b2 {name}: {out.runs[name]}")


PARTS = {"biology": part_biology, "donor_b2": part_donor_b2, "cases_b2": part_cases_b2}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.amendment5")
    ap.add_argument("--part", choices=[*PARTS, "all"], default="all")
    ap.add_argument("--only", default=None, metavar="REGEX", help="run only the runs whose tag matches (see amendment2)")
    ap.add_argument("--dataset", default=None, help="donor_b2: one data set only (one Kaggle session each)")
    args = ap.parse_args(argv)
    A2.SHARD = args.only
    for p in (PARTS if args.part == "all" else [args.part]):
        PARTS[p](args)


if __name__ == "__main__":
    main()
