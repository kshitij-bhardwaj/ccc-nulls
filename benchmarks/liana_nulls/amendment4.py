"""Amendment 4 (AMENDMENT_4.md): the leave-one-out t-test (LOO-t, PR #11's T2) through calibration, power, real
data and biology, and Null B2 with ligands shuffled.

    python -m benchmarks.liana_nulls.amendment4 [--part calibration|power|real|all]

LOO-t needs only each row's score, so the tools run with B = 10. Results: benchmarks/results/liana_nulls/
amendment4_*.json, each with AMENDMENT_4's SHA-256. Every part resumes where it stopped.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from ccc_nulls.nullb_type_restd import p_values

from .amendment3 import TISSUES, VALIDITY_REPS, shuffle_genes
from .common import AMENDMENT_DATASETS, CORE_DATASETS, OUT, amendment_hash, cpdb_database, log, rows_dir, save
from .core import ALPHA, KEY, add_q, bh_over, keyed
from .stages import _rows, kept_vs_dropped, mean_ci, spatial_truth

B_SCORE = 10            # permutations run only to obtain the tools' scores; LOO-t uses none
MIN_REF = 3
TOOLS = ("cellchat", "cellphonedb")


def save4(name: str, obj: dict):
    save(f"amendment4_{name}", {"amendment_4_sha256": amendment_hash(4), **obj})


def load4(name: str) -> dict:
    f = OUT / f"amendment4_{name}.json"
    return json.loads(f.read_text()) if f.exists() else {}


def loo_t(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """LOO-t p (1 where untestable) and BH calls over the tool's family."""
    fam = df.family.to_numpy(bool)
    p = p_values(df.assign(tested=fam), score="score", min_ref=MIN_REF)
    p = np.where(np.isfinite(p) & fam, p, 1.0)
    return p, bh_over(p, fam) <= ALPHA


def make_tool(name: str):
    from .tools import CellChat, CellPhoneDB
    return CellChat(B=B_SCORE, run_r=False) if name == "cellchat" else CellPhoneDB(cpdb_database(), B=B_SCORE, run_native=False)


def score_rows(tool: str, adata, groupby: str, tag: str) -> pd.DataFrame:
    df, _ = make_tool(tool).run(adata, groupby, rows_dir() / "amendment4" / tag / tool, log=lambda *a: None)
    df = add_q(keyed(df))
    for c in ("source", "target"):
        df[c] = df[c].astype(str)
    return df


def false_rows(rows: pd.DataFrame, side: str, shuffled: set) -> np.ndarray:
    col = "receptor_complex" if side == "receptor" else "ligand_complex"
    return np.array([all(g in shuffled for g in str(cx).split("_")) for cx in rows[col]])


def fdp(calls: np.ndarray, false: np.ndarray) -> float:
    return float((calls & false).sum() / calls.sum()) if calls.any() else 0.0


# ------------------------------------------------------------------ calibration
def part_calibration(args):
    from .data import decoy_universe, load
    from .extra import b2_columns

    res = load4("calibration").get("results", {})
    for name in TISSUES:
        if all(t in res.get(name, {}) for t in TOOLS):
            continue
        adata, groupby, _ = load(name, extra_genes=decoy_universe(name))
        res.setdefault(name, {})
        for tool in TOOLS:
            if tool in res[name]:
                continue
            stored = pd.read_pickle(rows_dir() / f"{name}_{tool}.pkl")[KEY + ["family"]]
            for c in ("source", "target"):
                stored[c] = stored[c].astype(str)
            stored = stored[stored.ligand_complex.notna() & stored.receptor_complex.notna()].reset_index(drop=True)
            lig = {g for cx in stored.ligand_complex for g in str(cx).split("_")}
            rec = {g for cx in stored.receptor_complex for g in str(cx).split("_")}
            runs = {"receptor": [], "ligand": [], "labels": []}
            for rep in range(VALIDITY_REPS):
                for scen, genes, seed in (("receptor", rec - lig, 5000), ("ligand", lig - rec, 6000)):
                    a, shuffled = shuffle_genes(adata, genes, np.random.default_rng(seed + rep))
                    df = score_rows(tool, a, groupby, f"{name}_{scen}{rep}")
                    _, c_loo = loo_t(df)
                    b2 = b2_columns(a, groupby, stored)
                    runs[scen].append({"loo_t": dict(fdp=fdp(c_loo, false_rows(df, scen, shuffled)), calls=int(c_loo.sum())),
                                       "mdt": dict(fdp=fdp(b2.call_b2.to_numpy(bool), false_rows(stored, scen, shuffled)),
                                                   calls=int(b2.call_b2.sum()))})
                    del a
                rng = np.random.default_rng(4000 + rep)
                a = adata.copy()
                a.obs[groupby] = pd.Categorical(rng.permutation(a.obs[groupby].astype(str).to_numpy()),
                                                categories=adata.obs[groupby].cat.categories)
                df = score_rows(tool, a, groupby, f"{name}_labels{rep}")
                _, c_loo = loo_t(df)
                b2 = b2_columns(a, groupby, stored)
                runs["labels"].append({"loo_t": dict(any_call=bool(c_loo.any()), calls=int(c_loo.sum())),
                                       "mdt": dict(any_call=bool(b2.call_b2.any()), calls=int(b2.call_b2.sum()))})
                del a
                log(f"calibration {name} {tool} rep {rep}: " + "; ".join(
                    f"{s} LOO-t {runs[s][-1]['loo_t'].get('fdp', runs[s][-1]['loo_t'].get('any_call'))} "
                    f"MDT {runs[s][-1]['mdt'].get('fdp', runs[s][-1]['mdt'].get('any_call'))}" for s in runs))
            summary = {}
            for t in ("loo_t", "mdt"):
                summary[t] = {s: dict(fdr=mean_ci([r[t]["fdp"] for r in runs[s]]), calls=float(np.mean([r[t]["calls"] for r in runs[s]])))
                              for s in ("receptor", "ligand")}
                summary[t]["labels"] = dict(p_any_call=float(np.mean([r[t]["any_call"] for r in runs["labels"]])),
                                            calls=float(np.mean([r[t]["calls"] for r in runs["labels"]])))
            res[name][tool] = {"summary": summary, "reps": VALIDITY_REPS, "runs": runs}
            save4("calibration", {"alpha": ALPHA, "results": res})
        del adata


# ------------------------------------------------------------------ power
def part_power(args):
    from . import power
    from .data import load

    out = load4("power")
    runs = out.get("runs", [])
    done = {(r["rep"], r["tool"]) for r in runs}
    adata, groupby, _ = load("tabula_sapiens_mammary")
    cands = power.candidate_pairs(adata, cpdb_database())
    f = 1.5
    for rep in range(5):
        if all((rep, t) in done for t in TOOLS):
            continue
        bg = power.background(adata, groupby, rep)
        plan = power.design(bg, groupby, cands, rep)
        a = power.plant(bg, groupby, plan, f, rep)
        for tool in TOOLS:
            if (rep, tool) in done:
                continue
            df = score_rows(tool, a, groupby, f"power_rep{rep}")
            _, df["call_loo_t"] = loo_t(df)
            truth = power.label_truth(df, plan)
            planted = set(zip(plan.source, plan.target, plan.ligand, plan.receptor))
            is_planted = np.array([k in planted for k in zip(df.source, df.target, df.ligand_complex, df.receptor_complex)])
            m = power.metrics(df, truth, is_planted, len(plan), {"loo_t": "call_loo_t"})
            runs.append(dict(rep=rep, effect=f, tool=tool, n_planted=len(plan), **m))
            log(f"power rep {rep} {tool}: LOO-t power {m['loo_t']['power']:.2f} fdp {m['loo_t']['fdp']:.3f}")
            summary = {t: dict(power=mean_ci([r["loo_t"]["power"] for r in runs if r["tool"] == t]),
                               fdr=mean_ci([r["loo_t"]["fdp"] for r in runs if r["tool"] == t]),
                               calls=float(np.mean([r["loo_t"]["calls"] for r in runs if r["tool"] == t])))
                       for t in TOOLS if any(r["tool"] == t for r in runs)}
            save4("power", {"effect": f, "summary": summary, "runs": runs})


# ------------------------------------------------------------------ real data and biology
def part_real(args):
    from . import hypothalamus
    out = {}
    for name in CORE_DATASETS + AMENDMENT_DATASETS:
        out[name] = {}
        if name == hypothalamus.NAME:
            truth = pd.read_csv(OUT / f"truth_{name}_spatial.csv")
        elif name in ("tnbc_wu2021", "allen_brain_atlas"):
            truth = spatial_truth(name)
        else:
            truth = None
        for tool in TOOLS:
            rows = _rows(name, tool)
            if rows is None:
                continue
            _, c = loo_t(rows)
            a = rows.call_null_a.to_numpy(bool)
            r = dict(calls=int(c.sum()), calls_null_a=int(a.sum()), calls_null_b=int(rows.call_null_b.sum()),
                     outside_null_a=int((c & ~a).sum()), jaccard_null_a=float((c & a).sum() / max((c | a).sum(), 1)))
            if truth is not None:
                r["primary_kept_vs_dropped"] = kept_vs_dropped(rows.assign(call_null_b=c), truth)
            out[name][tool] = r
            log(f"real {name} {tool}: LOO-t calls {r['calls']} (Null A {r['calls_null_a']}, Null B {r['calls_null_b']})"
                + (f"; OR {r['primary_kept_vs_dropped']['odds_ratio']:.2f}" if truth is not None else ""))
    save4("real", {"alpha": ALPHA, "results": out})


PARTS = {"calibration": part_calibration, "power": part_power, "real": part_real}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_nulls.amendment4")
    ap.add_argument("--part", choices=[*PARTS, "all"], default="all")
    args = ap.parse_args(argv)
    for p in (PARTS if args.part == "all" else [args.part]):
        PARTS[p](args)


if __name__ == "__main__":
    main()
