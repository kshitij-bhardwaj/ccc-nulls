"""AMENDMENT_2.md: held-out confirmation (A2.1), harder power (A2.2), donor-aware permutations (A2.3),
filter baselines (A2.4, inside every F1/F2 run), sensitivity (A2.5), p-value calibration (A2.6), a second
species (A2.7), biology multiplicity and detectable effects (A2.8), case studies (A2.9), run time (A2.10).

    python -m benchmarks.liana_final.amendment2 --part all
    python -m benchmarks.liana_final.amendment2 --part confirm

Every part writes benchmarks/results/liana_final/amendment2/<part>.json after each run and skips the
runs already there, so a stopped part resumes where it stopped.
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

from .evaluate import (DATASETS, SPATIAL, auroc, consensus_candidates, false_rows, load, side_only_genes, spatial_truth,
                       table, top)
from .pipeline import (ALPHA, BASELINE, HERE, N_PERMS, OUT, SEED, VARIANTS, amendment_hash, consensus, liana_view,
                       plan_hash, resource, rows_dir, wide)

OUT2 = OUT / "amendment2"
FILTERS = ["liana6_prop20", "liana6_prop30", "liana6_de"]
MEAN = ["liana6@mean", "swap_null_b@mean", "extra_null_b@mean"]
ALL = VARIANTS + FILTERS + MEAN
NULL_B = ["swap_null_b", "extra_null_b"]
DONOR = {"tnbc_wu2021": "orig.ident", "tabula_sapiens_mammary": "donor_id", "tabula_sapiens_liver": "donor_id",
         "tabula_sapiens_lung": "donor_id"}
LUNG = "tabula_sapiens_lung"
MAMMARY = "tabula_sapiens_mammary"


def a2_dir():
    d = rows_dir() / "amendment2"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ------------------------------------------------------------------ results, resumable
class Results:
    """benchmarks/results/liana_final/amendment2/<part>.json: {"runs": {tag: record}, ...}."""

    def __init__(self, part: str):
        self.part, self.f = part, OUT2 / f"{part}.json"
        self.obj = json.loads(self.f.read_text()) if self.f.exists() else {"runs": {}}

    def __contains__(self, tag):
        return tag in self.obj["runs"]

    def add(self, tag: str, rec: dict):
        self.obj["runs"][tag] = rec
        self.write()

    def write(self, **extra):
        OUT2.mkdir(parents=True, exist_ok=True)
        self.obj.update(extra)
        head = {"plan_sha256": plan_hash(), "amendment_1_sha256": amendment_hash(1),
                "amendment_2_sha256": hashlib.sha256((HERE / "AMENDMENT_2.md").read_bytes()).hexdigest(),
                "written": time.strftime("%Y-%m-%d %H:%M")}
        self.f.write_text(json.dumps({**head, **{k: v for k, v in self.obj.items() if k not in head}}, indent=1, default=str))

    @property
    def runs(self):
        return self.obj["runs"]


# ------------------------------------------------------------------ one run: every ranking of A2
def type_stats(a, groupby: str, de: bool = True):
    """(detection share, up-regulated) per cell type and gene: types x genes. Up: one-vs-rest t-test
    (scanpy ``t-test_overestim_var``), log fold change > 0 and BH-adjusted p < 0.05."""
    from scipy import sparse
    lab = a.obs[groupby].astype(str).to_numpy()
    types = sorted(set(lab))
    X = sparse.csr_matrix(a.X)
    genes = a.var_names.astype(str)
    prop = pd.DataFrame(np.stack([np.asarray((X[lab == t] > 0).mean(0)).ravel() for t in types]), index=types, columns=genes)
    if not de:
        return prop, None
    import scanpy as sc
    b = a.copy()
    b.obs["_g"] = pd.Categorical(lab)
    sc.tl.rank_genes_groups(b, "_g", method="t-test_overestim_var", use_raw=False)
    up = pd.DataFrame(0.0, index=types, columns=genes)
    for t in types:
        d = sc.get.rank_genes_groups_df(b, t)
        up.loc[t, d.loc[(d.logfoldchanges > 0) & (d.pvals_adj < 0.05), "names"].astype(str).to_numpy()] = 1.0
    return prop, up


def per_row(tbl: pd.DataFrame, M: pd.DataFrame, cx_col: str, type_col: str) -> np.ndarray:
    """Each row's value of M (types x genes) for its complex in its type; the lowest subunit counts."""
    cxs = pd.Index(tbl[cx_col].astype(str).unique())
    C = pd.DataFrame({cx: M.reindex(columns=cx.split("_")).fillna(0.0).min(axis=1) for cx in cxs})
    return C.to_numpy()[C.index.get_indexer(tbl[type_col].astype(str)), C.columns.get_indexer(tbl[cx_col].astype(str))]


def add_filters(tbl: pd.DataFrame, a, groupby: str) -> pd.DataFrame:
    """A2.4: liana6 with rows failing a detection or DE filter moved below every other row."""
    prop, up = type_stats(a, groupby)
    lp, rp = per_row(tbl, prop, "ligand_complex", "source"), per_row(tbl, prop, "receptor_complex", "target")
    lu, ru = per_row(tbl, up, "ligand_complex", "source"), per_row(tbl, up, "receptor_complex", "target")
    base = tbl[f"specificity_rank__{BASELINE}"].to_numpy(float)
    keep = {"liana6_prop20": (lp >= 0.2) & (rp >= 0.2), "liana6_prop30": (lp >= 0.3) & (rp >= 0.3),
            "liana6_de": (lu >= 1) & (ru >= 1)}
    for v, ok in keep.items():
        tbl[f"specificity_rank__{v}"] = base + (~ok)
        tbl[f"call__{v}"] = tbl[f"specificity_rank__{v}"].to_numpy() <= ALPHA
    return tbl


def run_all(adata, groupby: str, species: str, extra=(), *, n_perms: int = N_PERMS, seed: int = SEED,
            expr_prop: float = 0.1, strata=None, filters: bool = True, b2: bool = True):
    """One LIANA+ run: every variant (RRA), three of them also with LIANA+'s ``mean`` aggregation, and the
    A2.4 filters. Returns (wide table, captured tools, seconds)."""
    from benchmarks.liana_nulls.extra import excluded_genes
    t0 = time.time()
    a = liana_view(adata, extra)
    lrs, store = lc.run(a.copy(), groupby, resource_name=resource(species), n_perms=n_perms, seed=seed,
                        expr_prop=expr_prop, strata=strata)
    if b2:
        lc.add_b2(store, adata, groupby, excluded_genes())
    V = lc.variants(lrs, store, nulls=lc.NULLS if b2 else ("null_a", "null_b"))
    Vm = lc.variants(lrs, store, nulls=("null_b",), aggregate_method="mean")
    V |= {f"{v}@mean": Vm[v] for v in ("liana6", "swap_null_b", "extra_null_b")}
    tbl = wide(V)
    if filters:
        tbl = add_filters(tbl, a, groupby)
    return tbl, store, round(time.time() - t0, 1)


def variants_in(tbl):
    return [v for v in ALL if f"specificity_rank__{v}" in tbl]


# ------------------------------------------------------------------ calibration (F1 as amended) + A2.6
def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def p_calibration(store: dict, scen: str, shuffled: set) -> dict:
    """A2.6, per tool: the p-values of the false rows (shuffled side) and BH at 5%."""
    out = {}
    for tool, d in store.items():
        f = false_rows(d, scen, shuffled)
        r = {"rows": int(len(d)), "false_rows": int(f.sum()), "false_share_all": float(f.mean())}
        for n in lc.NULLS:
            p = d[f"p_{n}"].to_numpy(float)
            ok = np.isfinite(p)
            pf = p[f & ok]
            q = np.full(len(p), np.nan)
            q[ok] = bh(p[ok])
            disc = q <= 0.05
            r[n] = dict(at_most={str(a): float((pf <= a).mean()) for a in (0.01, 0.05, 0.1)},
                        quantiles=np.quantile(pf, np.linspace(0, 1, 101)).round(5).tolist() if len(pf) else [],
                        bh_discoveries=int(disc.sum()), bh_false_share=float(f[disc].mean()) if disc.any() else None,
                        rows_with_p=int(ok.sum()))
        out[tool] = r
    return out


def calibration_run(adata, groupby, sp, extra, scen, rep, only, *, pcal=False, **kw) -> dict:
    from benchmarks.liana_nulls.amendment3 import shuffle_genes
    a, shuffled = shuffle_genes(adata, only[scen], np.random.default_rng((7000 if scen == "receptor" else 8000) + rep))
    tbl, store, secs = run_all(a, groupby, sp, extra, **kw)
    f = false_rows(tbl, scen, shuffled)
    r = {"scenario": scen, "rep": rep, "rows": int(len(tbl)), "false_share_all": float(f.mean()), "seconds": secs,
         "false_share_top": {v: float(f[top(tbl, v)].mean()) for v in variants_in(tbl)}}
    if pcal:
        r["p_calibration"] = p_calibration(store, scen, shuffled)
    return r


def calibration_part(part: str, name: str, *, reps=3, scenarios=("receptor", "ligand"), adata=None, groupby=None,
                     tag="", pcal=False, **kw):
    res = Results(part)
    if adata is None:
        adata, groupby, sp, extra = load(name)
    else:
        from benchmarks.liana_nulls.data import SPECIES
        sp, extra = SPECIES[name], set()
        if name.startswith("tabula_sapiens"):
            from benchmarks.liana_nulls.data import decoy_universe
            extra = decoy_universe(name)
    only = side_only_genes(sp)
    for rep in range(reps):
        for scen in scenarios:
            t = f"{name}{tag}|{scen}|{rep}"
            if t in res:
                continue
            r = calibration_run(adata, groupby, sp, extra, scen, rep, only, pcal=pcal, **kw)
            res.add(t, r | {"dataset": name, "setting": tag or "default"})
            log(f"{part} {t}: all {r['false_share_all']:.3f}; top 1% " + ", ".join(
                f"{v} {r['false_share_top'][v]:.3f}" for v in ("liana6", "swap_null_b", "extra_null_b", "liana6_de")
                if v in r["false_share_top"]))
    return res


# ------------------------------------------------------------------ power (F2 design and A2.2)
def background(adata, groupby: str, rep: int, n: int):
    """power.background with n cells (n = 10,000 is power.background exactly)."""
    rng = np.random.default_rng(1000 + rep)
    keep = np.sort(rng.choice(adata.n_obs, min(n, adata.n_obs), replace=False))
    bg = adata[keep].copy()
    bg.obs[groupby] = pd.Categorical(rng.permutation(bg.obs[groupby].astype(str).to_numpy()))
    vc = bg.obs[groupby].value_counts()
    return bg[bg.obs[groupby].isin(vc[vc >= 20].index)].copy()


def complex_candidates(adata, min_frac: float = 0.02, max_sub: int = 3) -> pd.DataFrame:
    """Pairs of LIANA+'s consensus resource with a complex (<= max_sub subunits a side), every subunit
    detected in >= min_frac of cells."""
    import liana as li
    from scipy import sparse
    res = li.rs.select_resource("consensus")
    res = res[res.ligand.str.contains("_") | res.receptor.str.contains("_")][["ligand", "receptor"]].drop_duplicates()
    X = sparse.csc_matrix(adata.X)
    frac = pd.Series(np.diff(X.indptr) / X.shape[0], index=adata.var_names.astype(str))
    ok = []
    for lig, rec in zip(res.ligand, res.receptor):
        L, R = lig.split("_"), rec.split("_")
        ok.append(len(L) <= max_sub and len(R) <= max_sub and not set(L) & set(R)
                  and all(frac.get(g, 0) >= min_frac for g in L + R))
    return res[np.array(ok)].reset_index(drop=True)


def hard_design(bg, groupby: str, singles: pd.DataFrame, complexes: pd.DataFrame, rep: int) -> pd.DataFrame:
    """A2.2: 30 single, 10 unequal and 10 complex pairs, disjoint genes, random (s, t), s != t."""
    from benchmarks.liana_nulls.power import pick_pairs
    rng = np.random.default_rng(2000 + rep)
    s = pick_pairs(singles, 40, rng).assign(cls=["single"] * 30 + ["unequal"] * 10)
    used = set(s.ligand) | set(s.receptor)
    rows = []
    for i in rng.permutation(len(complexes)):
        lig, rec = complexes.ligand.iat[i], complexes.receptor.iat[i]
        g = set(lig.split("_")) | set(rec.split("_"))
        if g & used:
            continue
        used |= g
        rows.append((lig, rec))
        if len(rows) == 10:
            break
    plan = pd.concat([s, pd.DataFrame(rows, columns=["ligand", "receptor"]).assign(cls="complex")], ignore_index=True)
    types = np.array(sorted(bg.obs[groupby].astype(str).unique()))
    st = np.array([rng.choice(types, 2, replace=False) for _ in range(len(plan))])
    return plan.assign(source=st[:, 0], target=st[:, 1])


def plant(bg, groupby: str, plan: pd.DataFrame, f: float, rep: int):
    """power.plant for every subunit; an unequal pair's receptor is raised by 1 + (f - 1) / 2. For single
    pairs this is power.plant exactly (same order, same random numbers)."""
    from scipy import sparse

    from benchmarks.liana_nulls.power import raise_gene
    rng = np.random.default_rng(3000 + rep)
    a = bg.copy()
    X = sparse.csc_matrix(a.X, dtype=np.float32)
    var = pd.Index(a.var_names.astype(str))
    labels = a.obs[groupby].astype(str).to_numpy()
    cls = plan.cls if "cls" in plan else pd.Series("single", index=plan.index)
    for p, c in zip(plan.itertuples(), cls):
        for g in p.ligand.split("_"):
            raise_gene(X, var.get_loc(g), np.where(labels == p.source)[0], f, rng)
        for g in p.receptor.split("_"):
            raise_gene(X, var.get_loc(g), np.where(labels == p.target)[0], 1 + (f - 1) / 2 if c == "unequal" else f, rng)
    a.X = sparse.csr_matrix(X)
    return a


def truth(tbl: pd.DataFrame, plan: pd.DataFrame) -> np.ndarray:
    """power.label_truth with complexes: a side is raised when every subunit is raised in its type."""
    raised = {}
    for p in plan.itertuples():
        raised |= {g: p.source for g in p.ligand.split("_")} | {g: p.target for g in p.receptor.split("_")}

    def side(cx, ty):
        hit = np.array([[raised.get(g) == t for g in str(c).split("_")] for c, t in zip(cx, ty)], dtype=object)
        return np.array(["yes" if all(h) else ("partial" if any(h) else "no") for h in hit])
    L = side(tbl.ligand_complex.to_numpy(), tbl.source.astype(str).to_numpy())
    R = side(tbl.receptor_complex.to_numpy(), tbl.target.astype(str).to_numpy())
    lab = np.full(len(tbl), "null", dtype=object)
    lab[(L == "yes") ^ (R == "yes")] = "one_sided"
    lab[(L == "partial") | (R == "partial")] = "ambiguous"
    lab[(L == "yes") & (R == "yes")] = "true"
    return lab


def power_measures(tbl: pd.DataFrame, plan: pd.DataFrame) -> dict:
    lab = truth(tbl, plan)
    key = list(zip(tbl.source.astype(str), tbl.target.astype(str), tbl.ligand_complex, tbl.receptor_complex))
    cls = plan.cls if "cls" in plan else pd.Series("single", index=plan.index)
    which = {(p.source, p.target, p.ligand, p.receptor): c for p, c in zip(plan.itertuples(), cls)}
    planted = np.array([which.get(k) for k in key], dtype=object)
    true, one = lab == "true", lab == "one_sided"
    out = {"planted_rows_present": {c: int((planted == c).sum()) for c in sorted(set(cls))},
           "n_true_rows": int(true.sum()), "n_one_sided_rows": int(one.sum()), "rows": int(len(tbl))}
    k = int(true.sum())
    for v in variants_in(tbl):
        t = top(tbl, v)
        order = np.argsort(tbl[f"specificity_rank__{v}"].to_numpy(), kind="stable")[:k]
        out[v] = {"precision_at_k": float(true[order].mean()) if k else float("nan")}
        for c in sorted(set(cls)):
            m = planted == c
            out[v][c] = dict(power_top=float(t[m].sum() / (cls == c).sum()), auroc=auroc(tbl, v, m, one))
    return out


def power_part(part: str, name: str, *, effects, reps, design: str = "f2", n_cells: int = 10_000, tag: str = "",
               rep_list=None, **kw):
    from benchmarks.liana_nulls import power
    res = Results(part)
    adata, groupby, sp, extra = load(name)
    base = adata[:, ~adata.var_names.isin(list(extra))]
    singles = consensus_candidates(base)
    complexes = complex_candidates(base) if design == "hard" else None
    for rep in (rep_list if rep_list is not None else range(reps)):
        tags = [f"{name}{tag}|{design}|f{f}|{rep}" for f in effects]
        if all(t in res for t in tags):
            continue
        bg = background(adata, groupby, rep, n_cells)
        plan = hard_design(bg, groupby, singles, complexes, rep) if design == "hard" else power.design(bg, groupby, singles, rep)
        for f, t in zip(effects, tags):
            if t in res:
                continue
            a = plant(bg, groupby, plan, f, rep)
            tbl, store, secs = run_all(a, groupby, sp, extra, **kw)
            r = power_measures(tbl, plan) | {"effect": f, "rep": rep, "seconds": secs, "dataset": name,
                                             "setting": tag or "default", "design": design}
            res.add(t, r)
            log(f"{part} {t}: " + ", ".join(f"{v} {r[v]['single']['power_top']:.2f}/{r[v]['single']['auroc']:.3f}"
                                            for v in ("liana6", "swap_null_b", "extra_null_b", "liana6_de")))
            del a, tbl, store
    return res


# ------------------------------------------------------------------ parts
def part_confirm(args):
    """A2.1: the frozen rule on Tabula Sapiens Lung."""
    calibration_part("confirm_calibration", LUNG, reps=3)
    power_part("confirm_power", LUNG, effects=(1.5, 2.0), reps=5)


def part_hard_power(args):
    """A2.2 (with A2.4's filters in every run)."""
    power_part("hard_power", MAMMARY, effects=(1.1, 1.25, 1.5), reps=5, design="hard")


def within_donor_labels(adata, groupby: str, donor: str, rng):
    a = adata.copy()
    lab = a.obs[groupby].astype(str).to_numpy().copy()
    d = a.obs[donor].astype(str).to_numpy()
    for k in np.unique(d):
        i = np.where(d == k)[0]
        lab[i] = lab[i][rng.permutation(len(i))]
    a.obs[groupby] = pd.Categorical(lab)
    return a


def part_donor(args):
    """A2.3 (a) labels shuffled within donor: tool-level type-I error under global and within-donor
    permutations; (b) real labels: swap_null_b from within-donor permutations."""
    from benchmarks.liana_nulls.data import SPECIES
    res = Results("donor_null")
    for name, donor in DONOR.items():
        adata, groupby, sp, extra = load(name)
        a0 = liana_view(adata, extra)
        for rep in range(3):
            for scheme in ("global", "donor"):
                t = f"{name}|{scheme}|{rep}"
                if t in res:
                    continue
                t0 = time.time()
                a = within_donor_labels(a0, groupby, donor, np.random.default_rng(9000 + rep))
                strata = a.obs[donor].astype(str) if scheme == "donor" else None
                _, store = lc.run(a, groupby, resource_name=resource(sp), n_perms=N_PERMS, seed=SEED, strata=strata,
                                  tools_only=True)
                r = {"dataset": name, "scheme": scheme, "rep": rep, "seconds": round(time.time() - t0, 1),
                     "n_donors": int(a0.obs[donor].nunique())}
                for tool, d in store.items():
                    d = d[(d.source != d.target) & (d.magnitude > 0)]
                    r[tool] = {"rows": int(len(d))} | {f"{n}_le_{al}": float((d[f"p_{n}"] <= al).mean())
                                                       for n in ("null_a", "null_b") for al in (0.01, 0.05)}
                res.add(t, r)
                log(f"donor_null {t}: " + ", ".join(f"{tool} A {r[tool]['null_a_le_0.05']:.3f} B {r[tool]['null_b_le_0.05']:.3f}"
                                                   for tool in lc.TOOLS))
        del adata, a0
    real = Results("donor_real")
    for name in ("tnbc_wu2021", "tabula_sapiens_mammary", "tabula_sapiens_liver"):
        if name in real:
            continue
        adata, groupby, sp, extra = load(name)
        tbl, store, secs = run_all(adata, groupby, sp, extra, strata=adata.obs[DONOR[name]].astype(str),
                                   filters=False, b2=False)
        tbl.to_pickle(a2_dir() / f"donor_{name}_wide.pkl")
        g = table(name)
        m = g.merge(tbl[lc.KEY + ["specificity_rank__swap_null_b", "call__swap_null_b"]], on=lc.KEY, suffixes=("", "_donor"))
        r = {"rows": int(len(m)), "seconds": secs}
        k = max(1, int(np.ceil(0.01 * len(m))))
        tg = set(np.argsort(m["specificity_rank__swap_null_b"].to_numpy(), kind="stable")[:k])
        td = set(np.argsort(m["specificity_rank__swap_null_b_donor"].to_numpy(), kind="stable")[:k])
        cg, cd = m["call__swap_null_b"].to_numpy(bool), m["call__swap_null_b_donor"].to_numpy(bool)
        r["swap_null_b"] = dict(top1_jaccard=len(tg & td) / len(tg | td), calls_global=int(cg.sum()), calls_donor=int(cd.sum()),
                                call_jaccard=float((cg & cd).sum() / (cg | cd).sum()))
        real.add(name, r)
        log(f"donor_real {name}: {r}")
        del adata, tbl, store


def part_recal(args):
    """A2.6: F1 on Mammary and Liver again, same seeds, keeping the tool-level p-values."""
    for name in ("tabula_sapiens_mammary", "tabula_sapiens_liver"):
        calibration_part("recalibration", name, reps=3, pcal=True)


def part_sensitivity(args):
    """A2.5 on Tabula Sapiens Mammary (permutations: the hypothalamus; see the README)."""
    from benchmarks.liana_nulls.data import load as load_dataset
    adata, groupby, sp, extra = load(MAMMARY)
    # seeds, real labels: against the saved run (seed 1337)
    res = Results("sensitivity_seed")
    base = table(MAMMARY)
    for seed in (1, 2):
        t = f"seed{seed}"
        if t in res:
            continue
        tbl, _, secs = run_all(adata, groupby, sp, extra, seed=seed, filters=False, b2=False)
        res.add(t, stability(base, tbl) | {"seconds": secs})
        log(f"sensitivity seed {seed}: {res.runs[t]}")
    # expr_prop: calibration rep 0 (both sides), power f = 1.5 reps 0-1
    for ep in (0.05, 0.2):
        calibration_part("sensitivity_calibration", MAMMARY, reps=1, adata=adata, groupby=groupby, tag=f"@expr_prop={ep}",
                         expr_prop=ep)
        power_part("sensitivity_power", MAMMARY, effects=(1.5,), reps=2, tag=f"@expr_prop={ep}", expr_prop=ep)
    # cells: 50% and 25% per type
    for frac in (0.5, 0.25):
        rng = np.random.default_rng(0)
        keep = np.concatenate([rng.choice(i, int(np.ceil(frac * len(i))), replace=False)
                               for i in (np.where(adata.obs[groupby].to_numpy() == t)[0] for t in adata.obs[groupby].cat.categories)])
        sub = adata[np.sort(keep)].copy()
        calibration_part("sensitivity_calibration", MAMMARY, reps=1, scenarios=("receptor",), adata=sub, groupby=groupby,
                         tag=f"@cells={frac}")
        power_part("sensitivity_power", MAMMARY, effects=(1.5,), reps=1, n_cells=int(10_000 * frac), tag=f"@cells={frac}")
        del sub
    del adata
    # resolution: coarse and fine labels
    from benchmarks.liana_nulls.data import clean, decoy_universe
    ex = decoy_universe(MAMMARY)
    for col in ("broad_cell_class", "free_annotation"):
        a, _, _ = load_dataset(MAMMARY, extra_genes=ex, obs_cols=(col,))
        vc = a.obs[col].value_counts()
        a = a[a.obs[col].isin(vc[vc >= 5].index)].copy()
        a.obs[col] = clean(a.obs[col])
        calibration_part("sensitivity_calibration", MAMMARY, reps=1, adata=a, groupby=col, tag=f"@labels={col}")
        del a
    # permutations, 1,000 vs 5,000: the hypothalamus (Mammary's 5,000 x 86,000-row permutation tensor
    # needs about 7 GB)
    H = "moffitt2018_hypothalamus"
    res = Results("sensitivity_permutations")
    if "real" not in res:
        adata, groupby, sp, extra = load(H)
        tbl, _, secs = run_all(adata, groupby, sp, extra, n_perms=5000, filters=False, b2=False)
        res.add("real", stability(table(H), tbl) | {"seconds": secs})
        del adata, tbl
    adata, groupby, sp, extra = load(H)
    for B in (1000, 5000):
        calibration_part("sensitivity_permutations_calibration", H, reps=1, scenarios=("receptor",), adata=adata,
                         groupby=groupby, tag=f"@n_perms={B}", n_perms=B, filters=False)


def stability(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    """Spearman of specificity_rank and Jaccard of the top 1%, rows in both tables."""
    from scipy.stats import spearmanr
    m = a.merge(b, on=lc.KEY, suffixes=("", "_2"))
    out = {"rows": int(len(m))}
    for v in ("liana6", "swap_null_b", "extra_null_b"):
        x, y = m[f"specificity_rank__{v}"].to_numpy(), m[f"specificity_rank__{v}_2"].to_numpy()
        k = max(1, int(np.ceil(0.01 * len(m))))
        tx, ty = set(np.argsort(x, kind="stable")[:k]), set(np.argsort(y, kind="stable")[:k])
        out[v] = dict(spearman=float(spearmanr(x, y).statistic), top1_jaccard=len(tx & ty) / len(tx | ty))
    return out


def part_mouse(args):
    """A2.7: F1 on the Allen brain atlas (mouse)."""
    calibration_part("mouse_calibration", "allen_brain_atlas", reps=3)


def part_timing(args):
    """A2.10: wall time of li.mt.rank_aggregate against this pipeline (six methods, Null B, B2)."""
    import liana as li
    res = Results("timing")
    for name in (MAMMARY, "moffitt2018_hypothalamus"):
        if name in res:
            continue
        adata, groupby, sp, extra = load(name)
        a = liana_view(adata, extra)
        t0 = time.time()
        li.mt.rank_aggregate(a.copy(), groupby, resource_name=resource(sp), n_perms=N_PERMS, seed=SEED, use_raw=False,
                             inplace=False)
        t1 = time.time()
        consensus(adata, groupby, sp, extra)
        t2 = time.time()
        res.add(name, dict(liana_rank_aggregate=round(t1 - t0, 1), this_pipeline=round(t2 - t1, 1),
                           ratio=round((t2 - t1) / (t1 - t0), 2), n_cells=int(adata.n_obs)))
        log(f"timing {name}: {res.runs[name]}")


# ------------------------------------------------------------------ A2.8 biology, A2.9 cases (from saved tables)
def boot_or(K: np.ndarray, n_boot: int = 2000, seed: int = 0) -> dict:
    """Odds ratio of K's 2x2 (rows: units; columns kept+, kept-, dropped+, dropped-) with the cluster
    bootstrap of PLAN.md, plus a two-sided bootstrap p and the detectable odds ratio."""
    def odds(k):
        a, b, c, d = k + (0.5 if (k == 0).any() else 0)
        return (a / b) / (c / d)
    tot = K.sum(0)
    if len(K) == 0 or tot[:2].sum() == 0 or tot[2:].sum() == 0:
        return dict(odds_ratio=None)
    rng = np.random.default_rng(seed)
    lb = np.log([odds(K[rng.integers(0, len(K), len(K))].sum(0)) for _ in range(n_boot)])
    se = float(lb.std(ddof=1))
    p = max(2 * min((lb <= 0).mean(), (lb >= 0).mean()), 1 / n_boot)
    return dict(odds_ratio=float(odds(tot)), ci95=np.exp(np.percentile(lb, [2.5, 97.5])).tolist(), p=float(min(p, 1.0)),
                se_log=se, detectable=[float(np.exp(-2.80 * se)), float(np.exp(2.80 * se))],
                kept_calls=int(tot[0] + tot[1]), dropped_calls=int(tot[2] + tot[3]))


def spatial_K(tbl: pd.DataFrame, truth_df: pd.DataFrame, v: str) -> np.ndarray:
    t = truth_df.set_index(["source", "target"]).colocalised
    r = tbl[(tbl.source != tbl.target) & tbl[f"call__{BASELINE}"].to_numpy(bool)]
    col = np.array([t.get(k, np.nan) for k in zip(r.source.astype(str), r.target.astype(str))], float)
    r = r.assign(coloc=col, kept=r[f"call__{v}"].to_numpy(bool))[np.isfinite(col)]
    g = r.groupby(["source", "target"]).agg(coloc=("coloc", "first"), kept=("kept", "sum"), n=("kept", "size"))
    g["dropped"] = g.n - g.kept
    return np.c_[g.kept * g.coloc, g.kept * (1 - g.coloc), g.dropped * g.coloc, g.dropped * (1 - g.coloc)].astype(float)


def downstream_K(r: pd.DataFrame, v: str) -> np.ndarray:
    kept = r[f"call__{v}"].to_numpy(bool)
    d = r.assign(k=kept.astype(int), dr=(~kept).astype(int), y=r["responds"].astype(int))
    g = d.groupby("unit")[["k", "dr", "y"]].apply(lambda x: pd.Series(dict(a=(x.k * x.y).sum(), b=(x.k * (1 - x.y)).sum(),
                                                                         c=(x.dr * x.y).sum(), e=(x.dr * (1 - x.y)).sum())))
    return g.to_numpy(float)


def part_biology(args):
    """A2.8: every F3 test with a bootstrap p, BH within variant and across all, and the detectable OR;
    swap_null_b_donor on TNBC (A2.3 b)."""
    others = [v for v in VARIANTS if v != BASELINE]
    tests = {}
    donor_f = a2_dir() / "donor_tnbc_wu2021_wide.pkl"
    donor = pd.read_pickle(donor_f)[lc.KEY + ["call__swap_null_b"]].rename(columns={"call__swap_null_b": "call__swap_null_b_donor"}) \
        if donor_f.exists() else None
    for name in SPATIAL:
        tbl, tr = table(name), spatial_truth(name)
        if donor is not None and name == "tnbc_wu2021":
            tbl = tbl.merge(donor, on=lc.KEY, how="left").fillna({"call__swap_null_b_donor": False})
        for v in others + (["swap_null_b_donor"] if f"call__swap_null_b_donor" in tbl else []):
            tests[(v, f"{name} spatial")] = boot_or(spatial_K(tbl, tr, v))
    for name in DATASETS:
        for res in ("cytosig", "progeny"):
            r = pd.read_pickle(rows_dir() / f"downstream_{name}_{res}.pkl")
            if donor is not None and name == "tnbc_wu2021":
                r = r.merge(donor, on=lc.KEY, how="left").fillna({"call__swap_null_b_donor": False})
            for v in others + (["swap_null_b_donor"] if "call__swap_null_b_donor" in r else []):
                tests[(v, f"{name} {res}")] = boot_or(downstream_K(r, v))
    ok = [k for k, t in tests.items() if t["odds_ratio"] is not None]
    q_all = bh(np.array([tests[k]["p"] for k in ok]))
    for k, q in zip(ok, q_all):
        tests[k]["q_all"] = float(q)
    for v in {k[0] for k in ok}:
        ks = [k for k in ok if k[0] == v]
        for k, q in zip(ks, bh(np.array([tests[k]["p"] for k in ks]))):
            tests[k]["q_variant"] = float(q)
    out = Results("biology")
    out.obj["tests"] = {f"{v} | {t}": r for (v, t), r in tests.items()}
    out.write()
    log(f"biology: {len(ok)} tests; q_all < 0.05: " + ", ".join(f"{v} {t}" for (v, t) in ok if tests[(v, t)]["q_all"] < 0.05))


def part_cases(args):
    """A2.9: limiting side and its breadth, dropped vs kept; the top dropped and promoted rows."""
    from scipy.stats import mannwhitneyu
    out = Results("cases")
    for name in DATASETS:
        if name in out:
            continue
        adata, groupby, sp, extra = load(name)
        prop, _ = type_stats(liana_view(adata, extra), groupby, de=False)
        del adata
        tbl = table(name)
        store = pd.read_pickle(rows_dir() / f"{name}_store.pkl")
        r = {}
        base = tbl[(tbl.source != tbl.target) & tbl[f"call__{BASELINE}"]]
        for tool, d in store.items():
            m = base.merge(d[lc.KEY + ["p_L_null_b", "p_R_null_b"]], on=lc.KEY)
            lim = np.where(m.p_L_null_b > m.p_R_null_b, "ligand", np.where(m.p_R_null_b > m.p_L_null_b, "receptor", "tie"))
            B = pd.DataFrame(np.repeat(prop.ge(0.1).mean(0).to_numpy()[None], len(prop), 0), index=prop.index,
                             columns=prop.columns)            # breadth of each gene, the same in every type
            breadth_l, breadth_r = per_row(m, B, "ligand_complex", "source"), per_row(m, B, "receptor_complex", "target")
            breadth = np.where(lim == "ligand", breadth_l, np.where(lim == "receptor", breadth_r, np.maximum(breadth_l, breadth_r)))
            drop = ~m["call__swap_null_b"].to_numpy(bool)
            mw = mannwhitneyu(breadth[drop], breadth[~drop], alternative="greater") if drop.any() and (~drop).any() else None
            r[tool] = dict(kept=int((~drop).sum()), dropped=int(drop.sum()),
                           limiting_side_dropped={s: int(((lim == s) & drop).sum()) for s in ("ligand", "receptor", "tie")},
                           limiting_side_kept={s: int(((lim == s) & ~drop).sum()) for s in ("ligand", "receptor", "tie")},
                           breadth_dropped_median=float(np.median(breadth[drop])) if drop.any() else None,
                           breadth_kept_median=float(np.median(breadth[~drop])) if (~drop).any() else None,
                           mann_whitney_p=float(mw.pvalue) if mw else None)
            m = m.assign(limiting_side=lim, breadth=breadth.round(3))
            cols = lc.KEY + ["specificity_rank__liana6", "specificity_rank__swap_null_b", "p_L_null_b", "p_R_null_b",
                             "limiting_side", "breadth"]
            m[drop].sort_values("specificity_rank__liana6").head(20)[cols].to_csv(OUT2 / f"cases_{name}_{tool}_dropped.csv", index=False)
        new = tbl[(tbl.source != tbl.target) & tbl["call__swap_null_b"] & ~tbl[f"call__{BASELINE}"]]
        new.sort_values("specificity_rank__swap_null_b").head(20)[lc.KEY + ["specificity_rank__liana6", "specificity_rank__swap_null_b"]] \
            .to_csv(OUT2 / f"cases_{name}_promoted.csv", index=False)
        r["promoted_rows"] = int(len(new))
        out.add(name, r)
        log(f"cases {name}: " + ", ".join(f"{t} dropped {x['dropped']} breadth {x['breadth_dropped_median']} vs {x['breadth_kept_median']} "
                                         f"p {x['mann_whitney_p']}" for t, x in r.items() if isinstance(x, dict)))


PARTS = {"confirm": part_confirm, "hard_power": part_hard_power, "donor": part_donor, "recal": part_recal,
         "sensitivity": part_sensitivity, "mouse": part_mouse, "timing": part_timing, "biology": part_biology,
         "cases": part_cases}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.amendment2")
    ap.add_argument("--part", choices=[*PARTS, "report", "all"], default="all")
    args = ap.parse_args(argv)
    for p in (PARTS if args.part == "all" else [args.part] if args.part != "report" else []):
        PARTS[p](args)
    if args.part in ("report", "all"):
        from .report2 import write
        write()


if __name__ == "__main__":
    main()
