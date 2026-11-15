"""Stages F0-F3 of PLAN.md, as amended by AMENDMENT_1.md (F1 and F2 judged on rankings), each a
``part_*`` function of the command line."""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from benchmarks.liana_nulls.common import OUT as NULLS_OUT
from benchmarks.liana_nulls.common import log

from .pipeline import ALPHA, BASELINE, VARIANTS, consensus, fidelity, rows_dir, save, summary, wide

DATASETS = ["tnbc_wu2021", "allen_brain_atlas", "tabula_sapiens_mammary", "tabula_sapiens_liver", "moffitt2018_hypothalamus"]
SPATIAL = ["tnbc_wu2021", "allen_brain_atlas", "moffitt2018_hypothalamus"]
TISSUES = ["tabula_sapiens_mammary", "tabula_sapiens_liver"]
CAL_REPS, POWER_REPS, EFFECTS = 3, 5, (1.5, 2.0)


def load(name: str):
    """(adata, groupby, species, decoy genes): Tabula Sapiens carries its decoy universe for B2."""
    from benchmarks.liana_nulls.data import SPECIES, decoy_universe
    from benchmarks.liana_nulls.data import load as load_dataset

    extra = decoy_universe(name) if name.startswith("tabula_sapiens") else set()   # non-database genes only
    adata, groupby, _ = load_dataset(name, extra_genes=extra)
    return adata, groupby, SPECIES[name], extra


def table(name: str) -> pd.DataFrame:
    return pd.read_pickle(rows_dir() / f"{name}_wide.pkl")


# ------------------------------------------------------------------ the five data sets (+ F0)
def part_run(args):
    for name in args.datasets:
        adata, groupby, sp, extra = load(name)
        log(f"run {name}: {adata.n_obs:,} cells, {adata.obs[groupby].nunique()} types, {adata.n_vars:,} genes")
        V, store, meta = consensus(adata, groupby, sp, extra)
        tbl = wide(V)
        tbl.to_pickle(rows_dir() / f"{name}_wide.pkl")
        pd.to_pickle(store, rows_dir() / f"{name}_store.pkl")
        f0 = fidelity(adata, groupby, sp, V, store, extra, standalone=name == "moffitt2018_hypothalamus")
        s = summary(tbl, store)
        log(f"  F0 passed: {f0['passed']}; calls " + ", ".join(f"{v} {s['calls'][v]['calls']:,}" for v in VARIANTS))
        save(f"run_{name}", {"dataset": name, "meta": meta, "fidelity": f0, **s})
        del adata, V, store, tbl


# ------------------------------------------------------------------ rankings (AMENDMENT_1.md)
TOP = 0.01


def top(tbl: pd.DataFrame, v: str, frac: float = TOP) -> np.ndarray:
    """The rows in the top ``frac`` of variant v's ranking (specificity_rank, ascending; ties by row order)."""
    k = max(1, int(np.ceil(frac * len(tbl))))
    out = np.zeros(len(tbl), bool)
    out[np.argsort(tbl[f"specificity_rank__{v}"].to_numpy(), kind="stable")[:k]] = True
    return out


def auroc(tbl: pd.DataFrame, v: str, pos: np.ndarray, neg: np.ndarray) -> float:
    """P(a positive row outranks a negative row) in variant v's ranking; ties count one half."""
    from scipy.stats import rankdata
    if not pos.any() or not neg.any():
        return float("nan")
    r = rankdata(-tbl[f"specificity_rank__{v}"].to_numpy()[pos | neg])       # high = ranked first
    lab = pos[pos | neg]
    n1, n0 = lab.sum(), (~lab).sum()
    return float((r[lab].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


# ------------------------------------------------------------------ F1 calibration
def false_rows(tbl: pd.DataFrame, side: str, shuffled: set) -> np.ndarray:
    col = "receptor_complex" if side == "receptor" else "ligand_complex"
    return np.array([all(g in shuffled for g in str(cx).split("_")) for cx in tbl[col]])


def side_only_genes(species: str = "human") -> dict:
    import liana as li
    res = li.rs.select_resource("consensus" if species == "human" else f"{species}consensus")
    lig = {g for x in res.ligand for g in str(x).split("_")}
    rec = {g for x in res.receptor for g in str(x).split("_")}
    return {"receptor": rec - lig, "ligand": lig - rec}


def part_calibration(args):
    """F1 as amended: false share in the top 1% of each ranking, receptor- or ligand-only genes shuffled."""
    from benchmarks.liana_nulls.amendment3 import shuffle_genes
    from benchmarks.liana_nulls.stages import mean_ci

    only = side_only_genes("human")
    out = {}
    for name in [t for t in TISSUES if t in args.datasets]:      # --dataset: one tissue per machine
        adata, groupby, sp, extra = load(name)
        runs = {"receptor": [], "ligand": []}
        for rep in range(args.cal_reps):
            for scen in runs:
                t0 = time.time()
                a, shuffled = shuffle_genes(adata, only[scen], np.random.default_rng((7000 if scen == "receptor" else 8000) + rep))
                V, store, _ = consensus(a, groupby, sp, extra)
                tbl = wide(V)
                f = false_rows(tbl, scen, shuffled)
                r = {"rows": int(len(tbl)), "false_rows": int(f.sum()), "false_share_all": float(f.mean()),
                     "seconds": round(time.time() - t0, 1)}
                for v in VARIANTS:
                    t, c = top(tbl, v), tbl[f"call__{v}"].to_numpy(bool)
                    r[v] = dict(false_share_top=float(f[t].mean()), calls=int(c.sum()),
                                fdp_calls=float((c & f).sum() / c.sum()) if c.any() else 0.0)
                runs[scen].append(r)
                log(f"calibration {name} {scen} rep {rep}: all {r['false_share_all']:.3f}; top 1% " + ", ".join(
                    f"{v} {r[v]['false_share_top']:.3f}" for v in VARIANTS))
                del a, V, store, tbl
            out[name] = {scen: {v: dict(false_share_top=mean_ci([x[v]["false_share_top"] for x in rs]),
                                        fdp_calls=mean_ci([x[v]["fdp_calls"] for x in rs]))
                                for v in VARIANTS} | {"false_share_all": mean_ci([x["false_share_all"] for x in rs]), "runs": rs}
                         for scen, rs in runs.items()}
            save("calibration", {"reps": args.cal_reps, "top": TOP, "n_side_only_genes": {k: len(v) for k, v in only.items()},
                                 "results": out})
        del adata


# ------------------------------------------------------------------ F2 power
def consensus_candidates(adata, min_frac: float = 0.02) -> pd.DataFrame:
    """Single-gene pairs of LIANA+'s consensus resource, both genes detected in >= min_frac of cells."""
    import liana as li
    from scipy import sparse
    res = li.rs.select_resource("consensus")
    pairs = res[~res.ligand.str.contains("_") & ~res.receptor.str.contains("_")][["ligand", "receptor"]].drop_duplicates()
    X = sparse.csc_matrix(adata.X)
    frac = pd.Series(np.diff(X.indptr) / X.shape[0], index=adata.var_names.astype(str))
    ok = (frac.reindex(pairs.ligand).fillna(0).to_numpy() >= min_frac) & (frac.reindex(pairs.receptor).fillna(0).to_numpy() >= min_frac)
    return pairs[ok].reset_index(drop=True)


def part_power(args):
    """F2 as amended: pairs from LIANA+'s resource; rank-based power, AUROC and precision."""
    from benchmarks.liana_nulls import power
    from benchmarks.liana_nulls.stages import mean_ci

    adata, groupby, sp, extra = load("tabula_sapiens_mammary")
    cands = consensus_candidates(adata[:, ~adata.var_names.isin(list(extra))])
    results = []
    for rep in range(args.power_reps):
        bg = power.background(adata, groupby, rep)
        plan = power.design(bg, groupby, cands, rep)
        for f in EFFECTS:
            t0 = time.time()
            a = power.plant(bg, groupby, plan, f, rep)
            V, store, _ = consensus(a, groupby, sp, extra)
            tbl = wide(V)
            truth = power.label_truth(tbl, plan)
            planted = set(zip(plan.source, plan.target, plan.ligand, plan.receptor))
            is_planted = np.array([k in planted for k in zip(tbl.source, tbl.target, tbl.ligand_complex, tbl.receptor_complex)])
            true, one = truth == "true", truth == "one_sided"
            r = dict(rep=rep, effect=f, n_planted=len(plan), planted_rows_present=int(is_planted.sum()), n_true_rows=int(true.sum()),
                     n_one_sided_rows=int(one.sum()), rows=int(len(tbl)), seconds=round(time.time() - t0, 1))
            for v in VARIANTS:
                t = top(tbl, v)
                k = int(true.sum())
                order = np.argsort(tbl[f"specificity_rank__{v}"].to_numpy(), kind="stable")[:k]
                r[v] = dict(power_top=float(t[is_planted].sum() / len(plan)), auroc_planted_vs_one_sided=auroc(tbl, v, is_planted, one),
                            precision_at_k=float(true[order].mean()) if k else float("nan"),
                            one_sided_in_top=int((t & one).sum()))
            results.append(r)
            log(f"power rep {rep} f {f}: planted present {r['planted_rows_present']}/{len(plan)}; " + ", ".join(
                f"{v} {r[v]['power_top']:.2f}/{r[v]['auroc_planted_vs_one_sided']:.3f}" for v in VARIANTS))
            del a, V, store, tbl
        summ = {f: {v: {m: mean_ci([r[v][m] for r in results if r["effect"] == f])
                        for m in ("power_top", "auroc_planted_vs_one_sided", "precision_at_k", "one_sided_in_top")}
                    for v in VARIANTS} | {"planted_rows_present": float(np.mean([r["planted_rows_present"] for r in results if r["effect"] == f]))}
                for f in EFFECTS if any(r["effect"] == f for r in results)}
        save("power", {"design": dict(n_cells=power.N_CELLS, n_pairs=power.N_PAIRS, effects=list(EFFECTS), reps=args.power_reps,
                                      n_candidate_pairs=len(cands), candidates="LIANA+ consensus resource, single-gene pairs",
                                      top=TOP), "summary": summ, "runs": results})


# ------------------------------------------------------------------ F3 biology
def spatial_truth(name: str) -> pd.DataFrame:
    from benchmarks.liana_nulls.stages import spatial_truth as truth
    if name == "moffitt2018_hypothalamus":
        return pd.read_csv(NULLS_OUT / f"truth_{name}_spatial.csv")
    return truth(name)


def part_biology(args):
    from benchmarks.liana_nulls import downstream
    from benchmarks.liana_nulls.data import SPECIES, to_human
    from benchmarks.liana_nulls.data import load as load_dataset
    from benchmarks.liana_nulls.stages import kept_vs_dropped, pair_level

    others = [v for v in VARIANTS if v != BASELINE]
    out = {"spatial": {}, "downstream": {}}
    for name in SPATIAL:
        tbl, truth = table(name), spatial_truth(name)
        d = {}
        for v in others:
            k = kept_vs_dropped(tbl.assign(call_null_a=tbl[f"call__{BASELINE}"], call_null_b=tbl[f"call__{v}"]), truth)
            if not (k["kept_calls"] and k["dropped_calls"]):
                k |= dict(odds_ratio=None, ci95=None)
            d[v] = k
        d["pair_level_auprc"] = pair_level(tbl, truth, {v: f"call__{v}" for v in VARIANTS}, ref=BASELINE)
        out["spatial"][name] = d
        log(f"spatial {name}: " + ", ".join(f"{v} {d[v]['odds_ratio'] if d[v]['odds_ratio'] is None else round(d[v]['odds_ratio'], 2)}"
                                           for v in others))
    extra = set(downstream.cytosig().index) | set(downstream.progeny().index)
    for name in DATASETS:
        sp = SPECIES[name]
        adata, groupby, _ = load_dataset(name, extra_genes=extra)
        if sp != "human":
            adata, _ = to_human(adata, sp)
        prof = downstream.contrast_profiles(adata, groupby)
        del adata
        rows = downstream.human_rows(table(name).assign(call_null_a=lambda t: t[f"call__{BASELINE}"]), sp)
        out["downstream"][name] = {}
        for res in ("cytosig", "progeny"):
            r = downstream.truth_rows(rows, prof, res)
            r.to_pickle(rows_dir() / f"downstream_{name}_{res}.pkl")
            out["downstream"][name][res] = {v: downstream.odds_ratio(r, r[f"call__{v}"].to_numpy(bool), ~r[f"call__{v}"].to_numpy(bool))
                                            for v in others} | {"lines": int(len(r))}
            log(f"downstream {name} {res}: {len(r):,} lines; " + ", ".join(
                f"{v} {o['odds_ratio'] if o['odds_ratio'] is None else round(o['odds_ratio'], 2)}"
                for v, o in out["downstream"][name][res].items() if v != "lines"))
    save("biology", out)


# ------------------------------------------------------------------ decision
def decision() -> dict:
    """The decision rule (AMENDMENT_1.md, replacing PLAN.md's) for the judged variants, from the saved
    results (None where a stage is missing)."""
    import json

    from .pipeline import JUDGED, OUT
    J = lambda n: json.loads((OUT / f"{n}.json").read_text()) if (OUT / f"{n}.json").exists() else None
    cal, pw, bio = J("calibration"), J("power"), J("biology")
    out = {}
    for v in JUDGED:
        c1 = None if cal is None or len(cal["results"]) < len(TISSUES) else all(
            cal["results"][t][s][v]["false_share_top"][0] < cal["results"][t][s][BASELINE]["false_share_top"][0]
            for t in cal["results"] for s in ("receptor", "ligand"))
        c2 = None
        if pw is not None and "2.0" in pw["summary"]:
            m, b = pw["summary"]["2.0"][v], pw["summary"]["2.0"][BASELINE]
            c2 = m["power_top"][0] >= 0.8 and m["auroc_planted_vs_one_sided"][0] >= b["auroc_planted_vs_one_sided"][0]
        c3 = None
        if bio is not None:
            cis = [bio["spatial"][d][v]["ci95"] for d in bio["spatial"] if bio["spatial"][d][v]["ci95"]]
            c3 = any(lo > 1 for lo, hi in cis) and not any(hi < 1 for lo, hi in cis)
        verdict = None if None in (c1, c2, c3) else ("improves LIANA+" if c1 and c2 and c3 else
                                                      "valid, not biologically distinguishable" if c1 and c2 else "not supported")
        out[v] = dict(calibration=c1, power=c2, biology=c3, verdict=verdict)
    return out
