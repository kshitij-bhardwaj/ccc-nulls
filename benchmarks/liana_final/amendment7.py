"""AMENDMENT_7.md: Connectome through donors, sensitivity, case studies and Benchmark 1's tool-level tests.

    python -m benchmarks.liana_final.amendment7 --part all
    python -m benchmarks.liana_final.amendment7 --part donor --only '^tabula_sapiens_lung\\|global\\|0$'

AMENDMENT_6's switch is on (Connectome captured, T2 in every run, ``swap3_<n>`` scored). Results:
benchmarks/results/liana_final/amendment7/. Every part resumes where it stopped.
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
from . import amendment6 as A6
from . import pipeline as P

OUT7 = P.OUT / "amendment7"
NULLS4, SWAP3, SHOWN = A6.NULLS4, A6.SWAP3, A6.SHOWN


def head() -> dict:
    return {**A6.head(), "amendment_7_sha256": hashlib.sha256((P.HERE / "AMENDMENT_7.md").read_bytes()).hexdigest()}


class Results(A2.Results):
    def __init__(self, part: str):
        self.part, self.f = part, OUT7 / f"{part}.json"
        self.obj = json.loads(self.f.read_text()) if self.f.exists() else {"runs": {}}

    def write(self, **extra):
        OUT7.mkdir(parents=True, exist_ok=True)
        self.obj.update(extra)
        h = head()
        self.f.write_text(json.dumps({**h, **{k: v for k, v in self.obj.items() if k not in h}}, indent=1, default=str))


def real_table(name: str) -> pd.DataFrame:
    """AMENDMENT_6's table of the real labels (every variant, swap3 included)."""
    return pd.read_pickle(A6.a6_dir() / f"{name}_wide.pkl")


def stability(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    from scipy.stats import spearmanr
    m = a.merge(b, on=lc.KEY, suffixes=("", "_2"))
    out = {"rows": int(len(m))}
    k = max(1, int(np.ceil(0.01 * len(m))))
    for v in SHOWN:
        x, y = m[f"specificity_rank__{v}"].to_numpy(), m[f"specificity_rank__{v}_2"].to_numpy()
        tx, ty = set(np.argsort(x, kind="stable")[:k]), set(np.argsort(y, kind="stable")[:k])
        out[v] = dict(spearman=float(spearmanr(x, y).statistic), top1_jaccard=len(tx & ty) / len(tx | ty))
    return out


def switch_on():
    A6.switch_on()
    A2.Results, A2.stability, A2.table = Results, stability, real_table


def tools_run(adata, groupby, sp, extra, *, b2=True, strata=None, n_perms=A2.N_PERMS):
    """The three captured tools (no consensus), with B2 and T2."""
    from benchmarks.liana_nulls.extra import excluded_genes
    _, store = lc.run(P.liana_view(adata, extra), groupby, resource_name=P.resource(sp), n_perms=n_perms, seed=A2.SEED,
                      strata=strata, tools_only=True)
    if b2:
        lc.add_b2(store, adata, groupby, excluded_genes())
    return store


TOOL_NULLS = ("null_a", "null_a_plus1", "null_b", "b2", "t2")


def bh_calls(d: pd.DataFrame, n: str, B: int = A2.N_PERMS) -> np.ndarray:
    """BH at 5% over the rows with a p. ``null_a_plus1``: LIANA+'s p = #/B as (1 + #)/(B + 1), Benchmark 1's
    convention (Null B's p already has it)."""
    p = d["p_null_a"].to_numpy(float) * B if n == "null_a_plus1" else d[f"p_{n}"].to_numpy(float)
    if n == "null_a_plus1":
        p = (p + 1) / (B + 1)
    ok = np.isfinite(p)
    q = np.ones(len(p))
    q[ok] = A2.bh(p[ok])
    return q <= 0.05


# ------------------------------------------------------------------ Benchmark 2 parts
def part_donor(args):
    res = Results("donor_null")
    for name, donor in A2.DONOR.items():
        if all(f"{name}|{s}|{rep}" in res for rep in range(3) for s in ("global", "donor")):
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        for rep in range(3):
            for scheme in ("global", "donor"):
                t = f"{name}|{scheme}|{rep}"
                if t in res:
                    continue
                t0 = time.time()
                a = adata.copy()
                lab, dn, rng = a.obs[groupby].astype(str).to_numpy().copy(), a.obs[donor].astype(str).to_numpy(), \
                    np.random.default_rng(9000 + rep)
                for k in np.unique(dn):                     # A2.within_donor_labels, without a second copy
                    i = np.where(dn == k)[0]
                    lab[i] = lab[i][rng.permutation(len(i))]
                a.obs[groupby] = pd.Categorical(lab)
                strata = a.obs[donor].astype(str) if scheme == "donor" else None
                store = tools_run(a, groupby, sp, extra, b2=scheme == "global", strata=strata)
                del a
                r = {"dataset": name, "scheme": scheme, "rep": rep, "seconds": round(time.time() - t0, 1),
                     "n_donors": int(adata.obs[donor].nunique())}
                nulls = NULLS4 if scheme == "global" else ("null_a", "null_b", "t2")
                for tool, d in store.items():
                    d = d[(d.source != d.target) & (d.magnitude > 0)]
                    r[tool] = {"rows": int(len(d))} | {f"{n}_le_{al}": float((np.nan_to_num(d[f"p_{n}"].to_numpy(float), nan=1.0) <= al).mean())
                                                       for n in nulls for al in (0.01, 0.05)}
                res.add(t, r)
                log(f"donor {t}: " + ", ".join(f"{tool} " + " ".join(f"{n} {r[tool][f'{n}_le_0.05']:.3f}" for n in nulls)
                                               for tool in store))
                del store
        del adata


def part_sensitivity(args):
    """A2.5, unchanged, with AMENDMENT_6's switch."""
    A2.part_sensitivity(args)


def part_cases(args):
    from scipy.stats import mannwhitneyu
    out = Results("cases")
    for name in A2.DATASETS:
        if name in out:
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        prop, _ = A2.type_stats(P.liana_view(adata, extra), groupby, de=False)
        del adata
        tbl = real_table(name)
        m = tbl[(tbl.source != tbl.target) & tbl[f"call__{P.BASELINE}"]].reset_index(drop=True)
        B = pd.DataFrame(np.repeat(prop.ge(0.1).mean(0).to_numpy()[None], len(prop), 0), index=prop.index, columns=prop.columns)
        breadth = np.maximum(A2.per_row(m, B, "ligand_complex", "source"), A2.per_row(m, B, "receptor_complex", "target"))
        med = lambda x: float(np.median(x)) if len(x) else None
        r = {}
        for n in ("null_b", "b2", "t2"):
            drop, drop2 = ~m[f"call__swap3_{n}"].to_numpy(bool), ~m[f"call__swap_{n}"].to_numpy(bool)
            mw = mannwhitneyu(breadth[drop], breadth[~drop], alternative="greater") if drop.any() and (~drop).any() else None
            r[f"swap3_{n}"] = dict(kept=int((~drop).sum()), dropped=int(drop.sum()), breadth_dropped_median=med(breadth[drop]),
                                   breadth_kept_median=med(breadth[~drop]), mann_whitney_p=float(mw.pvalue) if mw else None,
                                   dropped_by_both=int((drop & drop2).sum()), dropped_swap3_only=int((drop & ~drop2).sum()),
                                   dropped_swap_only=int((~drop & drop2).sum()))
        out.add(name, r)
        log(f"cases {name}: {r}")


# ------------------------------------------------------------------ Benchmark 1 measures, tool level
def part_tool_labels(args):
    res = Results("tool_labels")
    for name in ("tabula_sapiens_mammary", "tabula_sapiens_liver"):
        if all(f"{name}|labels|{rep}" in res for rep in range(3)):
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        for rep in range(3):
            t = f"{name}|labels|{rep}"
            if t in res:
                continue
            a = adata.copy()
            a.obs[groupby] = pd.Categorical(np.random.default_rng(4000 + rep).permutation(a.obs[groupby].astype(str).to_numpy()))
            store = tools_run(a, groupby, sp, extra)
            del a
            r = {"dataset": name, "rep": rep}
            for tool, d in store.items():
                r[tool] = {"rows": int(len(d))} | {n: dict(any_call=bool(bh_calls(d, n).any()), calls=int(bh_calls(d, n).sum()))
                                                   for n in TOOL_NULLS}
            res.add(t, r)
            log(f"tool_labels {t}: " + ", ".join(f"{tool} " + " ".join(f"{n} {r[tool][n]['calls']}" for n in TOOL_NULLS) for tool in store))
        del adata


def part_tool_power(args):
    """Benchmark 1's power at tool level: f = 1.5, 5 replicates, B = 1,000 and 5,000 (as Benchmark 1 for BPT)."""
    from benchmarks.liana_nulls import power
    res = Results("tool_power")
    f, Bs = 1.5, (1000, 5000)
    if all(f"{A2.MAMMARY}|B{B}|f{f}|{rep}" in res for rep in range(5) for B in Bs):
        return
    adata, groupby, sp, extra = A2.load_data(A2.MAMMARY)
    singles = A2.consensus_candidates(adata[:, ~adata.var_names.isin(list(extra))])
    for rep in range(5):
        for B in Bs:
            t = f"{A2.MAMMARY}|B{B}|f{f}|{rep}"
            if t in res:
                continue
            bg = A2.background(adata, groupby, rep, 10_000)
            plan = power.design(bg, groupby, singles, rep)
            a = A2.plant(bg, groupby, plan, f, rep)
            store = tools_run(a, groupby, sp, extra, n_perms=B)
            del a
            planted = set(zip(plan.source, plan.target, plan.ligand, plan.receptor))
            r = {"rep": rep, "effect": f, "n_perms": B, "n_planted": int(len(plan))}
            for tool, d in store.items():
                d = d.drop_duplicates(lc.KEY).reset_index(drop=True)
                d = d.assign(**{f"call_{n}": bh_calls(d, n, B) for n in TOOL_NULLS})
                is_p = np.array([k in planted for k in zip(d.source.astype(str), d.target.astype(str), d.ligand_complex,
                                                           d.receptor_complex)])
                r[tool] = power.metrics(d, power.label_truth(d, plan), is_p, len(plan), {n: f"call_{n}" for n in TOOL_NULLS})
            res.add(t, r)
            log(f"tool_power {t}: " + ", ".join(f"{tool} " + " ".join(f"{n} {r[tool][n]['power']:.2f}/{r[tool][n]['fdp']:.2f}"
                                                                       for n in TOOL_NULLS) for tool in store))
            del store


def part_tool_real(args):
    res = Results("tool_real")
    for name in A2.DATASETS:
        if name in res:
            continue
        adata, groupby, sp, extra = A2.load_data(name)
        store = tools_run(adata, groupby, sp, extra)
        del adata
        r = {}
        for tool, d in store.items():
            d = d.drop_duplicates(lc.KEY).reset_index(drop=True)
            calls = {n: bh_calls(d, n) for n in TOOL_NULLS}
            r[tool] = {"rows": int(len(d))} | {f"calls_{n}": int(c.sum()) for n, c in calls.items()}
            if name in A2.SPATIAL:
                tbl = d[lc.KEY].assign(**{f"call__{P.BASELINE}": calls["null_a_plus1"]}, **{f"call__{n}": calls[n] for n in TOOL_NULLS[2:]})
                tr = A2.spatial_truth(name)
                r[tool]["spatial"] = {n: A2.boot_or(A2.spatial_K(tbl, tr, n)) for n in TOOL_NULLS[2:]}
        res.add(name, r)
        log(f"tool_real {name}: " + ", ".join(f"{t} " + " ".join(f"{n} {r[t][f'calls_{n}']:,}" for n in TOOL_NULLS) for t in r))
        del store
    tests = [(k, t, n) for k, r in res.runs.items() for t in r if isinstance(r[t], dict) and "spatial" in r[t]
             for n in r[t]["spatial"] if r[t]["spatial"][n].get("odds_ratio") is not None]
    if tests:
        for (k, t, n), q in zip(tests, A2.bh(np.array([res.runs[k][t]["spatial"][n]["p"] for k, t, n in tests]))):
            res.runs[k][t]["spatial"][n]["q_all"] = float(q)
        res.write()


PARTS = {"donor": part_donor, "sensitivity": part_sensitivity, "cases": part_cases, "tool_labels": part_tool_labels,
         "tool_power": part_tool_power, "tool_real": part_tool_real}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.amendment7")
    ap.add_argument("--part", choices=[*PARTS, "all"], default="all")
    ap.add_argument("--only", default=None, metavar="REGEX", help="run only the runs whose tag matches (see amendment2)")
    args = ap.parse_args(argv)
    A2.SHARD = args.only
    switch_on()
    for p in (PARTS if args.part == "all" else [args.part]):
        PARTS[p](args)


if __name__ == "__main__":
    main()
