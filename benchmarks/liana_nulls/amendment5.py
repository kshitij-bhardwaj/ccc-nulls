"""Amendment 5 (AMENDMENT_5.md): NATMI through every stage, with Null A, Null B, B2 and LOO-t.

    python -m benchmarks.liana_nulls.amendment5                         # core, calibration, power, biology
    python -m benchmarks.liana_nulls.amendment5 --part core --dataset tnbc_wu2021
    python -m benchmarks.liana_nulls.amendment5 --part core --no-native  # column 1 skipped (no NATMI install)

Parts:
  core         columns 1-3 (native NATMI's score, its Python recomputation with Null A, Null B), B2 and LOO-t on
               the five data sets; rows saved as data/liana_nulls/<dataset>_natmi.pkl, joined onto the cached
               liana_res when it exists. LIANA+'s consensus is not run again.
  calibration  Tabula Sapiens Mammary and Liver: receptor-only or ligand-only lrc2p genes shuffled (20 replicates
               each), labels shuffled (10); every null on the same replicates. Needs core on both tissues (B2
               reads the core rows, as in Amendments 3 and 4).
  power        the plan's planted data (replicates 0-4, f = 1.5, 2, 3, 4), every null.
  biology      spatial (TNBC, Allen brain, hypothalamus), cytokine (TNBC), CytoSig and PROGENy (all five).

Needs, besides the benchmark's environment: for column 1, the NATMI clone and its own environment
(tools/natmi/requirements-native.txt: NATMI_DIR, NATMI_PYTHON). Results: benchmarks/results/liana_nulls/
amendment5_*.json, each with AMENDMENT_5's SHA-256. Every part resumes where it stopped.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from .amendment3 import TISSUES, clopper, shuffle_genes
from .amendment4 import false_rows, fdp, loo_t
from .common import AMENDMENT_DATASETS, CORE_DATASETS, OUT, amendment_hash, cpdb_database, log, rows_dir, save
from .core import ALPHA, KEY, add_q, bh_over, keyed, write_back

TOOL = "natmi"
NULLS = ("null_a", "null_b", "b2", "loo_t")
CAL_REPS = {"receptor": 20, "ligand": 20, "labels": 10}
CAL_SEEDS = {"labels": 1000, "receptor": 2000, "ligand": 3000}      # tools/cellchat/null_b_calibration.py's
SPATIAL = ("tnbc_wu2021", "allen_brain_atlas")


def save5(name: str, obj: dict):
    save(f"amendment5_{name}", {"amendment_5_sha256": amendment_hash(5), **obj})


def load5(name: str) -> dict:
    f = OUT / f"amendment5_{name}.json"
    return json.loads(f.read_text()) if f.exists() else {}


def natmi_genes() -> set:
    from tools.natmi import adapter
    return adapter.database_genes()


def b2_excluded() -> set:
    """B2's decoys stay outside every database of the benchmark and lrc2p."""
    from .extra import excluded_genes
    return excluded_genes() | natmi_genes()


def make_tool(species: str, B: int, native: bool = True, seed: int | None = None):
    from .tools import NATMI
    return NATMI(B=B, run_native=native, species=species, seed=seed)


def with_loo_t(df: pd.DataFrame) -> pd.DataFrame:
    """LOO-t (Amendment 4) on NATMI's score over its family: p, BH q and calls."""
    p, c = loo_t(df)
    df["p_loo_t"], df["call_loo_t"] = p, c
    df["q_loo_t"] = bh_over(p, df.family.to_numpy(bool))
    return df


def call_columns(df: pd.DataFrame) -> dict:
    return {n: f"call_{n}" for n in NULLS if f"call_{n}" in df}


def str_keys(df: pd.DataFrame) -> pd.DataFrame:
    for c in ("source", "target"):
        df[c] = df[c].astype(str)
    return df


# ------------------------------------------------------------------ core (T0 and the real-data calls)
def part_core(args):
    from .core import run_dataset
    from .data import decoy_universe, load

    for name in args.dataset or CORE_DATASETS + AMENDMENT_DATASETS:
        f = rows_dir() / f"{name}_{TOOL}.pkl"
        if f.exists() and (OUT / f"amendment5_core_{name}.json").exists() and not args.force:
            log(f"core {name}: done")
            continue
        log("core:", name)
        genes = natmi_genes()
        adata, groupby, extras = load(name, extra_genes=genes)
        decoys = decoy_universe(name)                         # B2's decoys for Tabula Sapiens, as Amendment 3
        extra_adata = load(name, extra_genes=decoys | genes)[0] if decoys else None
        sp = extras["species"]
        rows, summary = run_dataset(adata, groupby, [make_tool(sp, args.B, native=not args.no_native)], rows_dir() / name,
                                    liana_aggregate=False, extra_adata=extra_adata, log=log)
        df = with_loo_t(rows[TOOL])
        s = summary[TOOL]
        s["calls"].update(b2_bh=int(df.call_b2.sum()), loo_t_bh=int(df.call_loo_t.sum()))
        s["outside_null_a"] = {n: int((df[c].to_numpy(bool) & ~df.call_null_a.to_numpy(bool)).sum())
                               for n, c in call_columns(df).items() if n != "null_a"}
        fid = s.get("fidelity")
        s["fidelity_passed"] = None if fid is None else bool(fid.get("identical"))
        df.to_pickle(f)
        lf = rows_dir() / f"{name}_liana_res.pkl"
        if lf.exists():
            lr, s["write_back"] = write_back(pd.read_pickle(lf), TOOL, df)
            lr.to_pickle(lf)
        save5(f"core_{name}", {"dataset": name, "groupby": groupby, "B": args.B, "native_run": not args.no_native,
                               "n_cells": summary["n_cells"], "n_types": summary["n_types"], TOOL: s})
        log(f"  calls: {s['calls']}; fidelity: {fid}")
        del adata, extra_adata, rows, df


# ------------------------------------------------------------------ calibration (T1, Amendments 3-4)
def calibration_summary(runs: dict) -> dict:
    from .stages import mean_ci
    out = {}
    for n in NULLS:
        out[n] = {s: dict(fdr=mean_ci([r[n]["fdp"] for r in runs[s]]), calls=float(np.mean([r[n]["calls"] for r in runs[s]])),
                          reps=len(runs[s]))
                  for s in ("receptor", "ligand") if runs.get(s)}
        if runs.get("labels"):
            out[n]["labels"] = dict(p_any_call=clopper([r[n]["any_call"] for r in runs["labels"]]),
                                    calls=float(np.mean([r[n]["calls"] for r in runs["labels"]])), reps=len(runs["labels"]))
    return out


def part_calibration(args):
    from tools.natmi import adapter

    from .data import decoy_universe, load
    from .extra import b2_columns

    out = load5("calibration")
    res = out.get("results", {})
    lig_only, rec_only = adapter.side_genes()
    excluded = b2_excluded()
    for name in TISSUES:
        core_f = rows_dir() / f"{name}_{TOOL}.pkl"
        if not core_f.exists():
            log(f"calibration {name}: no core rows; run --part core --dataset {name} first")
            continue
        stored = str_keys(pd.read_pickle(core_f)[KEY + ["family"]])
        stored = stored[stored.ligand_complex.notna() & stored.receptor_complex.notna()].reset_index(drop=True)
        runs = res.setdefault(name, {}).setdefault("runs", {s: [] for s in CAL_REPS})
        if all(len(runs[s]) >= n for s, n in CAL_REPS.items()):
            continue
        adata, groupby, _ = load(name, extra_genes=decoy_universe(name) | natmi_genes())
        var = set(adata.var_names.astype(str))
        genes_for = {"receptor": rec_only & var, "ligand": lig_only & var}
        res[name].update(n_ligand_only_genes=len(genes_for["ligand"]), n_receptor_only_genes=len(genes_for["receptor"]))
        for scen, n in CAL_REPS.items():
            for rep in range(len(runs[scen]), n):
                rng = np.random.default_rng(CAL_SEEDS[scen] + rep)
                if scen == "labels":
                    a, shuffled = adata.copy(), None
                    a.obs[groupby] = pd.Categorical(rng.permutation(adata.obs[groupby].astype(str).to_numpy()))
                else:
                    a, shuffled = shuffle_genes(adata, genes_for[scen], rng)
                t0 = time.time()
                df, _ = make_tool("human", args.B, native=False, seed=rep + 1).run(
                    a, groupby, rows_dir() / "amendment5" / f"{name}_{scen}{rep}", log=lambda *x: None)
                df = with_loo_t(add_q(str_keys(keyed(df))))
                b2 = b2_columns(a, groupby, stored, excluded=excluded)
                calls = {c: df[f"call_{c}"].to_numpy(bool) for c in ("null_a", "null_b", "loo_t")}
                if shuffled is None:
                    r = {c: dict(any_call=bool(v.any()), calls=int(v.sum())) for c, v in calls.items()}
                    r["b2"] = dict(any_call=bool(b2.call_b2.any()), calls=int(b2.call_b2.sum()))
                else:
                    fr, fs = false_rows(df, scen, shuffled), false_rows(stored, scen, shuffled)
                    r = {c: dict(fdp=fdp(v, fr), calls=int(v.sum())) for c, v in calls.items()}
                    r["b2"] = dict(fdp=fdp(b2.call_b2.to_numpy(bool), fs), calls=int(b2.call_b2.sum()))
                r["seconds"] = round(time.time() - t0, 1)
                runs[scen].append(r)
                res[name]["summary"] = calibration_summary(runs)
                save5("calibration", {"alpha": ALPHA, "B": args.B, "reps": CAL_REPS, "seeds": CAL_SEEDS, "results": res})
                log(f"calibration {name} {scen} rep {rep}: " + "; ".join(
                    f"{c} {r[c].get('fdp', r[c].get('any_call'))} ({r[c]['calls']})" for c in NULLS))
                del a, df, b2
        del adata


# ------------------------------------------------------------------ power (T2)
def power_summary(results: list) -> dict:
    from .stages import mean_ci
    df = pd.DataFrame(results)
    out = {}
    for f, g in df.groupby("effect"):
        recs = g.to_dict("records")
        out[str(f)] = {m: dict(power=mean_ci([r[m]["power"] for r in recs]), fdr=mean_ci([r[m]["fdp"] for r in recs]),
                               one_sided_false_calls=float(np.mean([r[m]["one_sided_false_calls"] for r in recs])),
                               calls=float(np.mean([r[m]["calls"] for r in recs])))
                       for m in recs[0] if isinstance(recs[0][m], dict) and "power" in recs[0][m]}
        out[str(f)]["reps"] = len(recs)
        out[str(f)]["planted_in_lrc2p"] = float(np.mean([r["n_planted_in_lrc2p"] for r in recs]))
    return out


def part_power(args):
    from tools.natmi import adapter

    from . import power
    from .data import database_genes, decoy_universe, load
    from .extra import b2_columns
    from .tools import NATMI

    name = "tabula_sapiens_mammary"
    adata, groupby, _ = load(name, extra_genes=decoy_universe(name) | natmi_genes())
    tool_genes = adata.var_names.isin(database_genes("human"))           # the plan's candidates, unchanged
    cands = power.candidate_pairs(adata[:, tool_genes], cpdb_database())
    in_db = set(zip(*[adapter.load_pairs()[c] for c in ("ligand", "receptor")]))
    excluded = b2_excluded()
    out = load5("power")
    results = out.get("runs", [])
    done = {(r["rep"], r["effect"]) for r in results}
    log(f"power: {len(cands)} candidate pairs; {len(done)} runs done")
    for rep in range(args.power_reps):
        bg = plan = None
        for f in power.EFFECTS:
            if (rep, f) in done:
                continue
            if bg is None:
                bg = power.background(adata, groupby, rep)
                plan = power.design(bg, groupby, cands, rep)
            t0 = time.time()
            a = power.plant(bg, groupby, plan, f, rep)
            df, _ = make_tool("human", args.B, native=False).run(a, groupby, rows_dir() / "amendment5" / "power" / f"rep{rep}_f{f}")
            df = with_loo_t(add_q(str_keys(keyed(df))))
            for null in ("null_a", "null_b"):                                # Amendment 2's tool-rule columns
                df[f"call_{null}_tool_rule"] = NATMI.tool_rule(df[f"p_{null}"].fillna(1.0), df)
            b2 = b2_columns(a, groupby, df, excluded=excluded)
            for c in b2.columns:
                df[c] = b2[c].to_numpy()
            plan_in = plan[[(l, r) in in_db for l, r in zip(plan.ligand, plan.receptor)]]
            truth = power.label_truth(df, plan)
            planted = set(zip(plan_in.source, plan_in.target, plan_in.ligand, plan_in.receptor))
            is_planted = np.array([k in planted for k in zip(df.source, df.target, df.ligand_complex, df.receptor_complex)])
            cols = {"native_rule": "call_port_native", **call_columns(df),
                    "null_a_tool_rule": "call_null_a_tool_rule", "null_b_tool_rule": "call_null_b_tool_rule"}
            m = power.metrics(df, truth, is_planted, max(len(plan_in), 1), cols)
            results.append(dict(rep=rep, effect=f, tool=TOOL, n_planted=len(plan), n_planted_in_lrc2p=len(plan_in),
                                n_family=int(df.family.sum()), seconds=round(time.time() - t0, 1), **m))
            save5("power", {"design": dict(n_cells=power.N_CELLS, n_pairs=power.N_PAIRS, effects=list(power.EFFECTS),
                                           reps=args.power_reps, B=args.B, n_candidate_pairs=len(cands)),
                            "summary": power_summary(results), "runs": results})
            log(f"power rep {rep} f {f}: " + "; ".join(f"{n} {m[n]['power']:.2f}/{m[n]['fdp']:.3f}" for n in NULLS)
                + f" ({len(plan_in)} of {len(plan)} planted pairs in lrc2p)")
            del a, df, b2


# ------------------------------------------------------------------ biology (T3, Amendment 1)
def _kept_vs_dropped(rows, truth, col) -> dict:
    from .stages import kept_vs_dropped
    k = kept_vs_dropped(rows.assign(call_null_b=rows[col].to_numpy(bool)), truth)
    if not (k["kept_calls"] and k["dropped_calls"]):          # undefined; the 0.5 correction would print a number
        k |= dict(odds_ratio=None, ci95=None)
    return k


def part_biology(args):
    from . import downstream, hypothalamus
    from .data import SPECIES, load, to_human
    from .stages import _rows, openproblems_metric, pair_level, spatial_truth

    out = load5("biology").get("results", {})
    extra = set(downstream.cytosig().index) | set(downstream.progeny().index)
    for name in args.dataset or CORE_DATASETS + AMENDMENT_DATASETS:
        rows = _rows(name, TOOL)
        if rows is None:
            log(f"biology {name}: no core rows; run --part core first")
            continue
        d = {}
        if name == hypothalamus.NAME:
            f = OUT / f"truth_{name}_spatial.csv"
            truth = pd.read_csv(f) if f.exists() else None
        elif name in SPATIAL:
            truth = spatial_truth(name)
        else:
            truth = None
        cols = call_columns(rows)
        if truth is not None:
            d["spatial_kept_vs_dropped"] = {n: _kept_vs_dropped(rows, truth, c) for n, c in cols.items() if n != "null_a"}
            d["pair_level_auprc"] = pair_level(rows, truth, {"native_rule": "call_port_native", **cols}
                                               | ({"liana_builtin": "call_liana"} if "call_liana" in rows else {}))
        q_cols = {"native_rule": ("p_port_native", "p"), **{n: (f"q_{n}", "q") for n in cols}}
        if name in SPATIAL:
            d["openproblems_source_target"] = openproblems_metric(rows, truth, ["source", "target"], q_cols)
        if name == "tnbc_wu2021":
            cyt = pd.read_csv(OUT / "truth_tnbc_wu2021_cytokine.csv")
            d["cytokine_ligand_target"] = openproblems_metric(rows, cyt.rename(columns={"response": "colocalised"}),
                                                              ["ligand", "target"], q_cols)
        sp = SPECIES[name]
        adata, groupby, _ = load(name, extra_genes=extra)
        if sp != "human":
            adata, _ = to_human(adata, sp)
        prof = downstream.contrast_profiles(adata, groupby)
        del adata
        rh = downstream.human_rows(rows, sp)
        for res in ("cytosig", "progeny"):
            r = downstream.truth_rows(rh, prof, res)
            r.to_pickle(rows_dir() / f"downstream_{name}_{TOOL}_{res}.pkl")
            e = {"null_b": downstream.endpoints(r)} if len(r) else {}
            for n, c in cols.items():
                if n not in ("null_a", "null_b") and len(r):
                    kept = r[c].to_numpy(bool)
                    e[n] = {"D1": downstream.odds_ratio(r, kept, ~kept)}
            d[res] = e
        out[name] = d
        save5("biology", {"alpha": ALPHA, "results": out})
        log(f"biology {name}: " + "; ".join(
            f"{n} OR {v['odds_ratio'] if v['odds_ratio'] is None else round(v['odds_ratio'], 2)}"
            for n, v in d.get("spatial_kept_vs_dropped", {}).items()))


PARTS = {"core": part_core, "calibration": part_calibration, "power": part_power, "biology": part_biology}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_nulls.amendment5")
    ap.add_argument("--part", nargs="+", choices=[*PARTS, "all"], default=["all"])
    ap.add_argument("--dataset", nargs="+")
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--power-reps", type=int, default=5)
    ap.add_argument("--no-native", action="store_true", help="skip column 1 (native NATMI); fidelity is not checked")
    ap.add_argument("--force", action="store_true", help="core: run again even where results exist")
    args = ap.parse_args(argv)
    for p in (PARTS if "all" in args.part else args.part):
        PARTS[p](args)


if __name__ == "__main__":
    main()
