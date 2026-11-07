"""Amendment 3 (AMENDMENT_3.md): extend, power_a3, b2_validity and ablation.

Every variant of the ablation is a call column computed by :func:`variant_calls` from one row table,
so T2, the validity checks and the five data sets score the same definitions.
"""
from __future__ import annotations

import json
import time
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats as sstats

from ccc_nulls import competitive, composite
from ccc_nulls.stats import bh

from .common import AMENDMENT_DATASETS, CORE_DATASETS, OUT, amendment_hash, cpdb_database, log, rows_dir, save
from .core import ALPHA, KEY, add_q, bh_over, keyed, write_back

POWER_RUNS = [(1.1, 1000), (1.2, 1000), (1.3, 1000), (1.5, 1000), (1.2, 5000), (1.5, 5000)]
VALIDITY_REPS = 10
TISSUES = ["tabula_sapiens_mammary", "tabula_sapiens_liver"]


def save3(name: str, obj: dict):
    save(name, {"amendment_3_sha256": amendment_hash(3), **obj})


def by_over(p, family):
    """Benjamini-Yekutieli q-values over the family (valid under any dependence)."""
    q = np.ones(len(p))
    fam = np.asarray(family, bool)
    m = int(fam.sum())
    if m:
        q[fam] = np.minimum(1.0, bh(np.asarray(p, float)[fam]) * np.sum(1.0 / np.arange(1, m + 1)))
    return q


def ligand_receptor_p(df: pd.DataFrame, tool: str) -> tuple[np.ndarray, np.ndarray]:
    """Null B's one-side p-values as (ligand side, receptor side). CellChat's L side is the ligand;
    CellPhoneDB's is partner a, the first type of "x|y"."""
    pL, pR = df.p_L_null_b.to_numpy(float), df.p_R_null_b.to_numpy(float)
    if tool == "cellchat":
        return pL, pR
    a_is_ligand = df.source.astype(str).to_numpy() == np.array([p.split("|", 1)[0] for p in df.cpdb_pair.astype(str)])
    return np.where(a_is_ligand, pL, pR), np.where(a_is_ligand, pR, pL)


def tool_rule(tool: str, p, df) -> np.ndarray:
    from .tools import CellChat, CellPhoneDB
    return (CellChat if tool == "cellchat" else CellPhoneDB).tool_rule(np.nan_to_num(np.asarray(p, float), nan=1.0), df)


def variant_calls(df: pd.DataFrame, tool: str) -> dict:
    """Call columns of every variant present in ``df`` (see AMENDMENT_3.md, A3.4)."""
    fam = df.family.to_numpy(bool)
    q = lambda p: bh_over(np.nan_to_num(np.asarray(p, float), nan=1.0), fam) <= ALPHA
    p_lig, p_rec = ligand_receptor_p(df, tool)
    out = {"null_a": df.call_null_a.to_numpy(bool), "null_b": df.call_null_b.to_numpy(bool),
           "null_b|ligand_only": q(p_lig), "null_b|receptor_only": q(p_rec),
           "null_b|by": by_over(df.p_null_b.fillna(1.0), fam) <= ALPHA,
           "null_b|tool_rule": tool_rule(tool, df.p_null_b, df), "null_a|tool_rule": tool_rule(tool, df.p_null_a, df)}
    if tool == "cellphonedb" and "n_pos_L_null_b" in df:
        testable = composite.testable(df.score, df.n_pos_L_null_b, df.n_pos_R_null_b)
        out["null_b|testable"] = (bh_over(df.p_null_b.fillna(1.0).to_numpy(), testable) <= ALPHA) & fam
    if "call_null_b_pr9" in df:
        out["null_b|pr9"] = df.call_null_b_pr9.to_numpy(bool)
        out["null_b|pr9_testable"] = df.call_null_b_pr9_testable.to_numpy(bool)
        out["null_b|pr9_by"] = df.call_null_b_pr9_by.to_numpy(bool)
    for v in ("b2", "b2_uncorrected", "b2_unmatched", "b2_joint"):
        if f"p_{v}" in df:
            out[v.replace("b2_", "b2|")] = q(df[f"p_{v}"])
    return out


# ------------------------------------------------------------------ extend
def part_extend(args):
    """B2 (both tools) and #9's Null B (CellPhoneDB) on the cached core rows of every data set; the
    rows and liana_res are rewritten with the new columns, columns 1-3 untouched."""
    from .data import SPECIES, decoy_universe, load
    from .extra import add_extra

    for name in args.dataset or CORE_DATASETS + AMENDMENT_DATASETS:
        adata, groupby, _ = load(name, extra_genes=decoy_universe(name))
        sp, res = SPECIES[name], {}
        for tool in ("cellchat", "cellphonedb"):
            f = rows_dir() / f"{name}_{tool}.pkl"
            if not f.exists():
                log(f"extend: no core rows for {name} / {tool}")
                continue
            df = pd.read_pickle(f)
            for c in ("source", "target"):
                df[c] = df[c].astype(str)
            t0 = time.time()
            log(f"extend {name} {tool}: {len(df):,} rows")
            df, meta = add_extra(adata, groupby, SimpleNamespace(name=tool, species=sp), df, rows_dir() / name / tool, sp, log=log)
            meta["seconds"] = round(time.time() - t0, 1)
            df.to_pickle(f)
            lf = rows_dir() / f"{name}_liana_res.pkl"
            if lf.exists():
                lr, cov = write_back(pd.read_pickle(lf), tool, df)
                lr.to_pickle(lf)
                meta["write_back"] = cov
            res[tool] = meta
            del df
        save3(f"extend_{name}", {"dataset": name, **res})
        del adata


# ------------------------------------------------------------------ T2, small effects and B = 5,000
def part_power_a3(args):
    from . import power
    from .data import database_genes, decoy_universe, load
    from .extra import b2_columns, pr9_columns
    from .tools import CellChat, CellPhoneDB

    name = "tabula_sapiens_mammary"
    adata, groupby, _ = load(name, extra_genes=decoy_universe(name))
    tool_genes = adata.var_names.isin(database_genes("human"))         # what T2 loaded before (power.json)
    cands = power.candidate_pairs(adata[:, tool_genes], cpdb_database())
    out_f = OUT / "power_a3.json"
    results = json.loads(out_f.read_text())["runs"] if out_f.exists() and not getattr(args, "force", False) else []
    done = {(r["rep"], r["effect"], r["B"], r["tool"]) for r in results}
    log(f"power_a3: {len(cands)} candidate pairs; {len(done)} runs already done")
    for f, B in POWER_RUNS:
        for rep in range(args.power_reps):
            if all((rep, f, B, t) in done for t in ("cellchat", "cellphonedb")):
                continue
            bg = power.background(adata, groupby, rep)
            plan = power.design(bg, groupby, cands, rep)
            a_full = power.plant(bg, groupby, plan, f, rep)
            a = a_full[:, a_full.var_names.isin(database_genes("human"))].copy()
            for tool in (CellChat(B=B, run_r=False), CellPhoneDB(cpdb_database(), B=B, run_native=False)):
                if (rep, f, B, tool.name) in done:
                    continue
                t0 = time.time()
                wd = rows_dir() / "power_a3" / f"rep{rep}_f{f}_B{B}" / tool.name
                df, meta = tool.run(a, groupby, wd, log=log)
                df = add_q(keyed(df))
                b2 = b2_columns(a_full, groupby, df)
                for c in b2.columns:
                    df[c] = b2[c].to_numpy()
                if tool.name == "cellphonedb":
                    pr, _ = pr9_columns(a, groupby, df, wd / "pr9", B=B)
                    for c in pr.columns:
                        df[c] = pr[c].to_numpy()
                truth = power.label_truth(df, plan)
                planted = set(zip(plan.source, plan.target, plan.ligand, plan.receptor))
                is_planted = np.array([k in planted for k in zip(df.source, df.target, df.ligand_complex, df.receptor_complex)])
                calls = variant_calls(df, tool.name)
                for k, v in calls.items():
                    df[f"_v_{k}"] = v
                m = power.metrics(df, truth, is_planted, len(plan), {k: f"_v_{k}" for k in calls})
                results.append(dict(rep=rep, effect=f, B=B, tool=tool.name, n_planted=len(plan), n_family=int(df.family.sum()),
                                    seconds=round(time.time() - t0, 1), **m))
                log(f"  f {f} B {B} rep {rep} {tool.name}: " + "; ".join(
                    f"{k} {m[k]['power']:.2f}/{m[k]['fdp']:.2f}" for k in ("null_a", "null_b", "null_b|tool_rule", "b2")))
                save3("power_a3", {"design": dict(runs=POWER_RUNS, reps=args.power_reps, n_candidate_pairs=len(cands)),
                                   "summary": power_a3_summary(results), "check_against_power_json": check_power(results),
                                   "runs": results})
                del df


def power_a3_summary(results) -> dict:
    from .stages import mean_ci
    out = {}
    df = pd.DataFrame(results)
    for (tool, f, B), g in df.groupby(["tool", "effect", "B"]):
        recs = g.to_dict("records")
        variants = [k for k in recs[0] if isinstance(recs[0][k], dict) and "power" in recs[0][k]]
        out.setdefault(tool, {})[f"f={f:g}, B={B}"] = {
            v: dict(power=mean_ci([r[v]["power"] for r in recs]), fdr=mean_ci([r[v]["fdp"] for r in recs]),
                    one_sided_false_calls=float(np.mean([r[v]["one_sided_false_calls"] for r in recs])),
                    calls=float(np.mean([r[v]["calls"] for r in recs])))
            for v in variants} | {"reps": len(recs)}
    return out


def check_power(results) -> dict:
    """f = 1.5, B = 1,000 must reproduce power.json's Null A and Null B columns run by run."""
    f = OUT / "power.json"
    if not f.exists():
        return {"checked": False}
    old = {(r["rep"], r["tool"]): r for r in json.loads(f.read_text())["runs"] if r["effect"] == 1.5}
    same, n = 0, 0
    for r in results:
        if r["effect"] == 1.5 and r["B"] == 1000 and (r["rep"], r["tool"]) in old:
            o = old[(r["rep"], r["tool"])]
            n += 1
            same += all(r[k][m] == o[k][m] for k in ("null_a", "null_b") for m in ("calls", "power", "fdp", "one_sided_false_calls"))
    return {"checked": n, "identical": same}


# ------------------------------------------------------------------ B2 validity
def shuffle_genes(adata, genes, rng):
    """A copy with the columns of ``genes`` moved to random cells, one permutation for all of them."""
    a = adata.copy()
    X = sparse.csc_matrix(a.X)
    idx = np.where(a.var_names.isin(list(genes)))[0]
    perm = rng.permutation(a.n_obs)
    cols = [X[:, j] for j in range(X.shape[1])]
    for j in idx:
        cols[j] = cols[j][perm]
    a.X = sparse.hstack(cols, format="csr")
    return a, set(a.var_names[idx])


def validity_one(a, groupby, rows: pd.DataFrame, shuffled: set | None) -> dict:
    """Side-level type-I error (labels shuffled) or FDR on rows with a shuffled receptor side."""
    from .extra import b2_columns
    b2 = b2_columns(a, groupby, rows)
    sides = b2.attrs.get("sides")
    out = {}
    calls = {v: b2[f"call_{v}"].to_numpy(bool) for v in ("b2", "b2_uncorrected", "b2_unmatched", "b2_joint")}
    if shuffled is None:
        g = sides[sides.gate]
        out["side_type1"] = {"b2": float((g.p < ALPHA).mean()), "b2_uncorrected": float((g.p_uncorrected < ALPHA).mean()),
                             "b2_unmatched": float((g.p_unmatched < ALPHA).mean()), "n_sides": int(len(g))}
        fam = rows.family.to_numpy(bool) & b2.gate_b2.to_numpy(bool)
        out["row_type1_joint"] = float((b2.p_b2_joint.to_numpy()[fam] < ALPHA).mean()) if fam.any() else float("nan")
        out["any_call"] = {v: bool(c.any()) for v, c in calls.items()}
        out["calls"] = {v: int(c.sum()) for v, c in calls.items()}
    else:
        false = np.array([all(g in shuffled for g in str(cx).split("_")) for cx in rows.receptor_complex])
        out["fdp"] = {v: float((c & false).sum() / c.sum()) if c.any() else 0.0 for v, c in calls.items()}
        out["calls"] = {v: int(c.sum()) for v, c in calls.items()}
        out["false_rows"] = int(false.sum())
    return out


def part_b2_validity(args):
    from .data import decoy_universe, load
    from .stages import mean_ci

    res = {}
    for name in TISSUES:
        adata, groupby, _ = load(name, extra_genes=decoy_universe(name))
        res[name] = {}
        for tool in ("cellchat", "cellphonedb"):
            f = rows_dir() / f"{name}_{tool}.pkl"
            if not f.exists():
                continue
            rows = pd.read_pickle(f)[KEY + ["family"]]
            for c in ("source", "target"):
                rows[c] = rows[c].astype(str)
            rows = rows[rows.ligand_complex.notna() & rows.receptor_complex.notna()].reset_index(drop=True)
            lig = {g for cx in rows.ligand_complex for g in str(cx).split("_")}
            rec_only = {g for cx in rows.receptor_complex for g in str(cx).split("_")} - lig
            labels, receptor = [], []
            for rep in range(VALIDITY_REPS):
                rng = np.random.default_rng(4000 + rep)
                a = adata.copy()
                a.obs[groupby] = pd.Categorical(rng.permutation(a.obs[groupby].astype(str).to_numpy()),
                                                categories=adata.obs[groupby].cat.categories)
                labels.append(validity_one(a, groupby, rows, None))
                del a
                a, shuffled = shuffle_genes(adata, rec_only, np.random.default_rng(5000 + rep))
                receptor.append(validity_one(a, groupby, rows, shuffled))
                del a
                log(f"b2_validity {name} {tool} rep {rep}: type-I {labels[-1]['side_type1']}; fdp {receptor[-1]['fdp']}")
            vs = ("b2", "b2_uncorrected", "b2_unmatched")
            res[name][tool] = {
                "labels_shuffled": {v: dict(side_type1=mean_ci([r["side_type1"][v] for r in labels]),
                                            p_any_call=clopper([r["any_call"][v] for r in labels])) for v in vs}
                | {"b2_joint": dict(row_type1=mean_ci([r["row_type1_joint"] for r in labels]),
                                    p_any_call=clopper([r["any_call"]["b2_joint"] for r in labels]))},
                "receptor_shuffled": {v: dict(fdr=mean_ci([r["fdp"][v] for r in receptor]),
                                              calls=float(np.mean([r["calls"][v] for r in receptor])))
                                      for v in vs + ("b2_joint",)},
                "reps": VALIDITY_REPS, "n_receptor_only_genes": len(rec_only), "runs": {"labels": labels, "receptor": receptor}}
            res[name][tool]["calibrated"] = calibrated(res[name][tool])
        del adata
        save3("b2_validity", {"alpha": ALPHA, "results": res})


def clopper(hits) -> dict:
    k, n = int(np.sum(hits)), len(hits)
    lo = sstats.beta.ppf(0.025, k, n - k + 1) if k else 0.0
    hi = sstats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return dict(k=k, n=n, p=k / n if n else float("nan"), ci95=[float(lo), float(hi)])


def calibrated(r: dict) -> dict:
    """AMENDMENT_3's criterion: the side type-I CI reaches 0.05 or below, and the FDR upper bound is <= 0.05."""
    out = {}
    for v in ("b2", "b2_uncorrected", "b2_unmatched"):
        t1 = r["labels_shuffled"][v]["side_type1"][1][0] <= ALPHA
        fdr = r["receptor_shuffled"][v]["fdr"][1][1] <= ALPHA
        out[v] = dict(type1_ok=bool(t1), fdr_ok=bool(fdr), calibrated=bool(t1 and fdr))
    return out


# ------------------------------------------------------------------ ablation on the five data sets
def _downstream_join(r: pd.DataFrame, df: pd.DataFrame, calls: dict) -> dict:
    """Variant calls onto the downstream lines ``r`` (a subset of df's rows, complexes maybe renamed)."""
    idc = "interaction_name" if "interaction_name" in df else "native_interaction_id"
    key = ["source", "target", idc, "ligand", "receptor"]
    v = df[key].copy()
    for k, c in calls.items():
        v[f"_v_{k}"] = c
    v = v.drop_duplicates(key)
    j = r[key].astype({"source": str, "target": str}).merge(v, on=key, how="left")
    return {k: j[f"_v_{k}"].fillna(False).to_numpy(bool) for k in calls}


def part_ablation(args):
    from . import downstream
    from .stages import _rows, kept_vs_dropped, spatial_truth

    out = {}
    for name in CORE_DATASETS + AMENDMENT_DATASETS:
        out[name] = {}
        if name == "moffitt2018_hypothalamus":
            f = OUT / f"truth_{name}_spatial.csv"
            truth = pd.read_csv(f) if f.exists() else None
        elif name in ("tnbc_wu2021", "allen_brain_atlas"):
            truth = spatial_truth(name)
        else:
            truth = None
        for tool in ("cellchat", "cellphonedb"):
            df = _rows(name, tool)
            if df is None:
                continue
            calls = variant_calls(df, tool)
            full = calls["null_b"]
            d = {}
            for v, c in calls.items():
                u = (c | full).sum()
                d[v] = dict(calls=int(c.sum()), jaccard_with_null_b=float((c & full).sum() / u) if u else float("nan"),
                            outside_null_a=int((c & ~calls["null_a"]).sum()))
                if truth is not None and v != "null_a":
                    k = kept_vs_dropped(df.assign(call_null_b=c), truth)
                    if not (k["kept_calls"] and k["dropped_calls"]):  # undefined; the 0.5 correction would print a number
                        k |= dict(odds_ratio=None, ci95=None)
                    d[v]["spatial_kept_vs_dropped"] = {x: k[x] for x in ("odds_ratio", "ci95", "kept_calls", "dropped_calls")}
                for res in ("cytosig", "progeny"):
                    pf = rows_dir() / f"downstream_{name}_{tool}_{res}.pkl"
                    if v == "null_a" or not pf.exists():
                        continue
                    r = pd.read_pickle(pf)
                    if len(r) == 0:
                        continue
                    kept = _downstream_join(r, df, {v: c})[v]
                    o = downstream.odds_ratio(r, kept, ~kept)
                    d[v][f"D1_{res}"] = {x: o.get(x) for x in ("odds_ratio", "ci95", "kept_calls", "dropped_calls")}
            out[name][tool] = d
            log(f"ablation {name} {tool}: " + "; ".join(f"{v} {d[v]['calls']}" for v in d))
            del df
    save3("ablation", out)
