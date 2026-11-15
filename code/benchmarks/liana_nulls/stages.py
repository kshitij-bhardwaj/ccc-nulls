"""T1 calibration, T2 power, T3 biology and the report (plan: PLAN.md); Amendment 1 (AMENDMENT_1.md)."""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats

from .common import AMENDMENT_DATASETS, CORE_DATASETS, OUT, amendment_hash, cpdb_database, log, rows_dir, save
from .core import ALPHA, add_q, keyed
from .amendment3 import part_ablation, part_b2_validity, part_extend, part_power_a3  # noqa: F401  (Amendment 3)

ROOT = Path(__file__).resolve().parents[2]
TISSUES = ["mammary", "liver"]
CALIBRATION_FILES = {
    "cellchat": (ROOT / "tools/cellchat/results/null_b_calibration_summary.json",
                 [sys.executable, "-m", "tools.cellchat.null_b_calibration", "--tissue", "{tissue}"]),
    "cellphonedb": (ROOT / "benchmarks/results/null_b_cellphonedb/calibration_{tissue}.json",
                    [sys.executable, "-m", "benchmarks.null_b_cellphonedb", "--db", "{db}", "--part", "calibration",
                     "--tissue", "{tissue}"]),
}


def mean_ci(x):
    x = np.asarray(x, float)
    if len(x) < 2 or np.all(x == x[0]):
        return float(x.mean()) if len(x) else float("nan"), [float(x.mean()) if len(x) else float("nan")] * 2
    h = sstats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return float(x.mean()), [float(x.mean() - h), float(x.mean() + h)]


# ------------------------------------------------------------------ T1 calibration
def part_calibration(args):
    """Each tool's calibration script, run through this entry point (``--recompute``) or read from its
    saved output. Both were run with B = 1,000, 20 replicates per one-sided scenario, 10 shuffled-label
    replicates, BH over the tool's frozen family."""
    out = {}
    for tool, (path, cmd) in CALIBRATION_FILES.items():
        out[tool] = {}
        for tissue in TISSUES:
            f = Path(str(path).format(tissue=tissue))
            if getattr(args, "recompute", False) or not f.exists():
                log(f"calibration: running {tool} on {tissue}")
                subprocess.run([c.format(tissue=tissue, db=cpdb_database()) for c in cmd], cwd=ROOT, check=True)
            d = json.loads(f.read_text())
            d = d[tissue] if tool == "cellchat" else d["scenarios"]
            out[tool][tissue] = {}
            for scen in ("receptor", "ligand", "labels"):
                row = {}
                for null in ("null_a", "null_b"):
                    s = d[scen][null] if tool == "cellchat" else d[scen]["summary"][null]
                    s = s["0.05"]
                    if scen == "labels":
                        row[null] = dict(p_any_call=s.get("p_any_false_call", s.get("p_any_call")),
                                         ci95=s.get("p_any_false_call_ci95", s.get("ci95")))
                    else:
                        fdr = s["fdr"] if isinstance(s["fdr"], (int, float)) else s["fdr"]["mean"]
                        ci = s["fdr_ci95"] if "fdr_ci95" in s else s["fdr"]["ci95"]
                        row[null] = dict(fdr=fdr, ci95=ci)
                out[tool][tissue][scen] = row
            out[tool][tissue]["source"] = str(f.relative_to(ROOT))
    save("calibration", {"alpha": ALPHA, "results": out})


# ------------------------------------------------------------------ T2 power
def part_power(args):
    from . import power
    from .data import load
    from .tools import CellChat, CellPhoneDB

    adata, groupby, _ = load("tabula_sapiens_mammary")
    cands = power.candidate_pairs(adata, cpdb_database())
    log(f"power: {len(cands)} candidate pairs")
    results, reps = [], range(args.power_reps)
    for rep in reps:
        bg = power.background(adata, groupby, rep)
        plan = power.design(bg, groupby, cands, rep)
        for f in power.EFFECTS:
            a = power.plant(bg, groupby, plan, f, rep)
            for tool in (CellChat(B=args.B, run_r=False), CellPhoneDB(cpdb_database(), B=args.B, run_native=False)):
                t0 = time.time()
                df, meta = tool.run(a, groupby, rows_dir() / "power" / f"rep{rep}_f{f}" / tool.name, log=log)
                df = add_q(keyed(df))
                for null in ("null_a", "null_b"):              # Amendment 2: the tool's own rule on each null
                    df[f"call_{null}_tool_rule"] = tool.tool_rule(df[f"p_{null}"].fillna(1.0), df)
                truth = power.label_truth(df, plan)
                planted = set(zip(plan.source, plan.target, plan.ligand, plan.receptor))
                is_planted = np.array([k in planted for k in zip(df.source, df.target, df.ligand_complex, df.receptor_complex)])
                m = power.metrics(df, truth, is_planted, len(plan),
                                  {"native_rule": "call_port_native", "null_a": "call_null_a", "null_b": "call_null_b",
                                   "null_a_tool_rule": "call_null_a_tool_rule", "null_b_tool_rule": "call_null_b_tool_rule"})
                results.append(dict(rep=rep, effect=f, tool=tool.name, n_planted=len(plan), n_family=int(df.family.sum()),
                                    seconds=round(time.time() - t0, 1), **m))
                log(f"  rep {rep} f {f} {tool.name}: A power {m['null_a']['power']:.2f} fdp {m['null_a']['fdp']:.3f} | "
                    f"B power {m['null_b']['power']:.2f} fdp {m['null_b']['fdp']:.3f}")
            save("power", {"amendment_2_sha256": amendment_hash(2), "design": dict(n_cells=power.N_CELLS, n_pairs=power.N_PAIRS, effects=list(power.EFFECTS),
                                          reps=args.power_reps, B=args.B, n_candidate_pairs=len(cands)),
                           "summary": power_summary(results), "runs": results})


def power_summary(results):
    df = pd.DataFrame(results)
    out = {}
    for (tool, f), g in df.groupby(["tool", "effect"]):
        out.setdefault(tool, {})[str(f)] = {
            m: dict(power=mean_ci([r[m]["power"] for r in g.to_dict("records")]),
                    fdr=mean_ci([r[m]["fdp"] for r in g.to_dict("records")]),
                    one_sided_false_calls=float(np.mean([r[m]["one_sided_false_calls"] for r in g.to_dict("records")])),
                    calls=float(np.mean([r[m]["calls"] for r in g.to_dict("records")])))
            for m in ("native_rule", "null_a", "null_b", "null_a_tool_rule", "null_b_tool_rule")
            if m in g.iloc[0]} | {"reps": int(len(g))}
    return out


# ------------------------------------------------------------------ T3 biology
def spatial_truth(dataset: str) -> pd.DataFrame:
    f = OUT / f"truth_{dataset}_spatial.csv"
    if not f.exists():
        from . import spatial
        from .data import load
        adata, groupby, extras = load(dataset)
        if dataset == "tnbc_wu2021":
            truth, info = spatial.tnbc_truth(adata, groupby, log=log)
            (OUT / f"truth_{dataset}_spatial.json").write_text(json.dumps(info, indent=1))
        else:
            truth = spatial.brain_truth(adata)
        truth.to_csv(f, index=False)
        if dataset == "tnbc_wu2021":
            extras["openproblems_truth"].to_csv(OUT / "truth_tnbc_wu2021_cytokine.csv", index=False)
    return pd.read_csv(f)


def kept_vs_dropped(rows: pd.DataFrame, truth: pd.DataFrame, n_boot: int = 2000, seed: int = 0) -> dict:
    """Primary endpoint: colocalisation odds ratio of Null A's calls that Null B keeps vs drops, with a
    bootstrap over (s, t) pairs (calls of one pair share its truth)."""
    t = truth.set_index(["source", "target"]).colocalised
    r = rows[(rows.source != rows.target) & rows.call_null_a.to_numpy(bool)]
    col = np.array([t.get(k, np.nan) for k in zip(r.source.astype(str), r.target.astype(str))], float)
    r = r.assign(coloc=col)[np.isfinite(col)]
    r = r.assign(kept=r.call_null_b.to_numpy(bool))
    g = r.groupby(["source", "target"]).agg(coloc=("coloc", "first"), kept=("kept", "sum"), n=("kept", "size"))
    g["dropped"] = g.n - g.kept
    K = np.c_[g.kept * g.coloc, g.kept * (1 - g.coloc), g.dropped * g.coloc, g.dropped * (1 - g.coloc)]

    def odds(k):
        a, b, c, d = k + (0.5 if (k == 0).any() else 0)
        return (a / b) / (c / d)
    tot = K.sum(0)
    rng = np.random.default_rng(seed)
    boot = np.array([odds(K[rng.integers(0, len(K), len(K))].sum(0)) for _ in range(n_boot)])
    return dict(odds_ratio=float(odds(tot)), ci95=[float(x) for x in np.percentile(boot, [2.5, 97.5])],
                kept_calls=int(tot[0] + tot[1]), kept_colocalised=int(tot[0]),
                dropped_calls=int(tot[2] + tot[3]), dropped_colocalised=int(tot[2]),
                frac_colocalised_kept=float(tot[0] / max(tot[0] + tot[1], 1)),
                frac_colocalised_dropped=float(tot[2] / max(tot[2] + tot[3], 1)), n_pairs=int(len(K)))


def pair_level(rows: pd.DataFrame, truth: pd.DataFrame, methods: dict, ref: str = "null_b", n_boot: int = 1000) -> dict:
    """(s, t)-level AUPRC, each pair scored by its number of calls; paired bootstrap over truth pairs."""
    from tools.cellchat.match_null_a import bootstrap_metrics
    joined = {}
    for name, col in methods.items():
        n = rows.loc[rows[col].fillna(False).to_numpy(bool)].groupby(["source", "target"]).size()
        s = np.array([n.get(k, 0) for k in zip(truth.source, truth.target)], float)
        joined[name] = (truth.colocalised.to_numpy(int), s)
    return bootstrap_metrics(joined, ref=ref, n_boot=n_boot)


def openproblems_metric(rows: pd.DataFrame, truth: pd.DataFrame, keys, q_cols: dict, ref="null_b", n_boot=1000) -> dict:
    """Open Problems' own metric: the tool's score kept where called, max per key, against the truth."""
    from ccc_nulls.evaluation import openproblems as opm
    from tools.cellchat.match_null_a import bootstrap_metrics
    truth = truth.rename(columns={"colocalised": "response"}).copy()
    truth["response"] = truth.response.astype(int)
    joined = {}
    for name, (col, rule) in q_cols.items():
        p = rows[col].fillna(1.0).to_numpy()
        pred = opm.predictions(rows, rows.score.to_numpy(), keys, p=p if rule == "q" else np.where(p < ALPHA, 0.0, 1.0),
                               alpha=ALPHA)
        gt = opm._join(pred, truth, keys)
        joined[name] = (gt.response.to_numpy().astype(int), gt.score.to_numpy())
    rnd = truth[list(keys)].copy()
    rnd["score"] = np.random.default_rng(1).random(len(rnd))
    gt = opm._join(rnd, truth, keys)
    joined["random"] = (gt.response.to_numpy().astype(int), gt.score.to_numpy())
    return bootstrap_metrics(joined, ref=ref, n_boot=n_boot)


def part_biology(args):
    out = {}
    methods = {"native_rule": "call_port_native", "null_a": "call_null_a", "null_b": "call_null_b"}
    for dataset in ("tnbc_wu2021", "allen_brain_atlas"):
        truth = spatial_truth(dataset)
        out[dataset] = {"spatial_truth": dict(pairs=int(len(truth)), colocalised=int(truth.colocalised.sum()))}
        for tool in ("cellchat", "cellphonedb"):
            f = rows_dir() / f"{dataset}_{tool}.pkl"
            if not f.exists():
                log(f"biology: no core rows for {dataset} / {tool}; run --part core first")
                continue
            rows = pd.read_pickle(f)
            for c in ("source", "target"):
                rows[c] = rows[c].astype(str)
            d = {"primary_kept_vs_dropped": kept_vs_dropped(rows, truth),
                 "pair_level_auprc": pair_level(rows, truth, methods | ({"liana_builtin": "call_liana"} if "call_liana" in rows else {})),
                 "openproblems_source_target": openproblems_metric(
                     rows, truth, ["source", "target"],
                     {"native_rule": ("p_port_native", "p"), "null_a": ("q_null_a", "q"), "null_b": ("q_null_b", "q")})}
            if dataset == "tnbc_wu2021":
                cyt = pd.read_csv(OUT / "truth_tnbc_wu2021_cytokine.csv")
                d["cytokine_ligand_target"] = openproblems_metric(
                    rows, cyt.rename(columns={"response": "colocalised"}), ["ligand", "target"],
                    {"native_rule": ("p_port_native", "p"), "null_a": ("q_null_a", "q"), "null_b": ("q_null_b", "q")})
            out[dataset][tool] = d
            pk = d["primary_kept_vs_dropped"]
            log(f"biology {dataset} {tool}: kept vs dropped OR {pk['odds_ratio']:.2f} {np.round(pk['ci95'], 2).tolist()}")
    save("biology", out)


# ------------------------------------------------------------------ Amendment 1
def _rows(dataset: str, tool: str):
    f = rows_dir() / f"{dataset}_{tool}.pkl"
    if not f.exists():
        return None
    rows = pd.read_pickle(f)
    for c in ("source", "target"):
        rows[c] = rows[c].astype(str)
    return rows


def amendment_a1() -> dict:
    """A1: the hypothalamus MERFISH neighbourhood truth, with the primary endpoint and the (s, t) AUPRC."""
    from . import hypothalamus
    name = hypothalamus.NAME
    f = OUT / f"truth_{name}_spatial.csv"
    if not f.exists():
        truth, info = hypothalamus.neighbourhood_truth(log=log)
        truth.to_csv(f, index=False)
        (OUT / f"truth_{name}_spatial.json").write_text(json.dumps(info, indent=1))
    truth = pd.read_csv(f)
    out = {"spatial_truth": dict(pairs=int(len(truth)), colocalised=int(truth.colocalised.sum()))}
    methods = {"native_rule": "call_port_native", "null_a": "call_null_a", "null_b": "call_null_b"}
    for tool in ("cellchat", "cellphonedb"):
        rows = _rows(name, tool)
        if rows is None:
            log(f"amendment A1: no core rows for {name} / {tool}; run --part core --dataset {name} first")
            continue
        out[tool] = {"primary_kept_vs_dropped": kept_vs_dropped(rows, truth), "pair_level_auprc": pair_level(rows, truth, methods)}
        pk = out[tool]["primary_kept_vs_dropped"]
        log(f"A1 {name} {tool}: kept vs dropped OR {pk['odds_ratio']:.2f} {np.round(pk['ci95'], 2).tolist()}")
    return out


def amendment_a2() -> dict:
    """A2: downstream-response truths (CytoSig, PROGENy) on every dataset with core rows."""
    from . import downstream
    from .data import SPECIES, load, to_human
    extra = set(downstream.cytosig().index) | set(downstream.progeny().index)
    out = {}
    for dataset in CORE_DATASETS + AMENDMENT_DATASETS:
        tools = {t: r for t in ("cellchat", "cellphonedb") if (r := _rows(dataset, t)) is not None}
        if not tools:
            continue
        sp = SPECIES[dataset]
        adata, groupby, _ = load(dataset, extra_genes=extra)
        if sp != "human":
            adata, _ = to_human(adata, sp)
        prof = downstream.contrast_profiles(adata, groupby)
        del adata
        out[dataset] = {}
        for tool, rows in tools.items():
            rows = downstream.human_rows(rows, sp)
            out[dataset][tool] = {}
            for res in ("cytosig", "progeny"):
                r = downstream.truth_rows(rows, prof, res)
                r.to_pickle(rows_dir() / f"downstream_{dataset}_{tool}_{res}.pkl")
                e = downstream.endpoints(r)
                out[dataset][tool][res] = e
                d1 = e["D1"]
                log(f"A2 {dataset} {tool} {res}: {e['n_lines']:,} lines, {e['n_units']} units; D1 OR "
                    f"{d1['odds_ratio'] if d1['odds_ratio'] is None else round(d1['odds_ratio'], 2)} {d1['ci95']}")
    return out


def part_amendment(args):
    save("amendment_1", {"A1_hypothalamus_spatial": amendment_a1(), "A2_downstream": amendment_a2()}, amendment=True)


# ------------------------------------------------------------------ report
def part_report(args):
    from .report import write
    write()
