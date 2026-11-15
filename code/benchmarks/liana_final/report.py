"""results/liana_final/README.md from the JSON files."""
from __future__ import annotations

import json

from .evaluate import DATASETS, SPATIAL, decision
from .pipeline import BASELINE, JUDGED, OUT, VARIANTS, amendment_hash, plan_hash

LABEL = {"liana5": "LIANA+ as shipped (5 methods)", "liana6": "**+ CellChat (baseline)**", "swap_null_b": "**swap: Null B**",
         "swap_b2": "swap: B2", "extra_null_a": "7th/8th: Null A", "extra_null_b": "**7th/8th: Null B**", "extra_b2": "7th/8th: B2"}


def _load(n):
    f = OUT / f"{n}.json"
    return json.loads(f.read_text()) if f.exists() else None


def _ci(x, d=3):
    return "" if not x else f"[{x[0]:.{d}f}, {x[1]:.{d}f}]"


def _or(o):
    if not o or o.get("odds_ratio") is None:
        return "–"
    lo, hi = o["ci95"]
    mark = "↑" if lo > 1 else ("↓" if hi < 1 else "")
    small = f" (k={o['kept_calls']})" if o.get("kept_calls", 20) < 20 else ""
    return f"{o['odds_ratio']:.2f} [{lo:.2f}, {hi:.2f}]{mark}{small}"


def write():
    L = ["# Our nulls inside LIANA+'s consensus: results", "",
         f"Plan: `benchmarks/liana_final/PLAN.md` (SHA-256 `{plan_hash()[:16]}…`), amended by `AMENDMENT_1.md` "
         f"(SHA-256 `{(amendment_hash(1) or '')[:16]}…`): F1 and F2 judged on the ranking, after one power replicate showed the "
         "call rule could not separate the variants. Run: `python -m benchmarks.liana_final`.",
         "Baseline: LIANA+'s five methods plus CellChat (liana6). A call (F0, F3) is specificity_rank ≤ 0.05.", ""]
    runs = {d: _load(f"run_{d}") for d in DATASETS}
    if any(runs.values()):
        L += ["## F0 fidelity and calls", "",
              "| Dataset | Rows | liana5 = rank_aggregate | " + " | ".join(LABEL[v] for v in VARIANTS) + " |",
              "|---|---|---|" + "---|" * len(VARIANTS)]
        for d, r in runs.items():
            if r:
                ok = "yes" if r["fidelity"]["passed"] else "**no**"
                cells = [f"{r['calls'][v]['calls']:,} ({r['calls'][v]['jaccard_with_liana6']:.2f})" for v in VARIANTS]
                L.append(f"| {d} | {r['rows']:,} | {ok} | " + " | ".join(cells) + " |")
        L += ["", "Calls (Jaccard with the baseline's calls).", ""]
    cal = _load("calibration")
    if cal:
        L += ["## F1 calibration (Tabula Sapiens, known nulls; AMENDMENT_1.md)", "",
              f"Share of rows with a shuffled side in the top {100 * cal['top']:.0f}% of each ranking, mean of {cal['reps']} "
              "replicates [95% t-CI]. The first line is the share among all rows. Lower is better.", "",
              "| Variant | " + " | ".join(f"{t.split('_')[-1]}, {s}s shuffled" for t in cal["results"] for s in ("receptor", "ligand")) + " |",
              "|---|" + "---|" * (2 * len(cal["results"]))]
        L.append("| all rows | " + " | ".join(f"{cal['results'][t][s]['false_share_all'][0]:.3f}"
                                              for t in cal["results"] for s in ("receptor", "ligand")) + " |")
        for v in VARIANTS:
            L.append(f"| {LABEL[v]} | " + " | ".join(f"{cal['results'][t][s][v]['false_share_top'][0]:.3f} {_ci(cal['results'][t][s][v]['false_share_top'][1])}"
                                                     for t in cal["results"] for s in ("receptor", "ligand")) + " |")
        L.append("")
    pw = _load("power")
    if pw:
        L += ["## F2 power (planted edges, Tabula Sapiens Mammary; AMENDMENT_1.md)", "",
              f"Pairs from LIANA+'s consensus resource. Power in the top {100 * pw['design']['top']:.0f}% / AUROC planted vs "
              f"one-sided / precision at k (k = true rows), mean of the replicates run "
              f"(planted rows present: " + ", ".join(f"f = {f}: {m['planted_rows_present']:.1f} of {pw['design']['n_pairs']}"
                                                      for f, m in pw["summary"].items()) + ").", "",
              "| Variant | " + " | ".join(f"f = {f}" for f in pw["summary"]) + " |", "|---|" + "---|" * len(pw["summary"])]
        for v in VARIANTS:
            L.append(f"| {LABEL[v]} | " + " | ".join(
                f"{m[v]['power_top'][0]:.2f} / {m[v]['auroc_planted_vs_one_sided'][0]:.3f} / {m[v]['precision_at_k'][0]:.2f}"
                for m in pw["summary"].values()) + " |")
        L.append("")
    bio = _load("biology")
    if bio:
        others = [v for v in VARIANTS if v != BASELINE]
        L += ["## F3 biology", "", "Odds ratio among the baseline's calls (s ≠ t), kept vs dropped by the variant "
              "(95% bootstrap CI; ↑ above 1, ↓ below 1; – where it keeps or drops none; k where it keeps fewer than 20).", "",
              "| Variant | " + " | ".join(f"{d} spatial" for d in bio["spatial"]) + " | "
              + " | ".join(f"{d} {r}" for d in bio["downstream"] for r in ("cytosig", "progeny")) + " |",
              "|---|" + "---|" * (len(bio["spatial"]) + 2 * len(bio["downstream"]))]
        for v in others:
            L.append(f"| {LABEL[v]} | " + " | ".join(_or(bio["spatial"][d][v]) for d in bio["spatial"]) + " | "
                     + " | ".join(_or(bio["downstream"][d][r][v]) for d in bio["downstream"] for r in ("cytosig", "progeny")) + " |")
        L += ["", "Pair-level AUPRC (pairs scored by their number of calls; difference from the baseline, paired bootstrap):", "",
              "| Dataset | " + " | ".join(LABEL[v] for v in VARIANTS) + " |", "|---|" + "---|" * len(VARIANTS)]
        for d in bio["spatial"]:
            m = bio["spatial"][d]["pair_level_auprc"]["methods"]
            L.append(f"| {d} | " + " | ".join(f"{m[v]['auprc']:.3f}" + (f" ({-m[v]['auprc_diff_ref_minus_this']:+.3f})"
                                                                         if "auprc_diff_ref_minus_this" in m[v] else "")
                                              for v in VARIANTS) + " |")
        L.append("")
    dec = decision()
    L += ["## Decision (AMENDMENT_1.md's rule)", "", "| Variant | Calibration | Power | Biology | Verdict |", "|---|---|---|---|---|"]
    for v in JUDGED:
        d = dec[v]
        fmt = lambda x: "pending" if x is None else ("yes" if x else "no")
        L.append(f"| {LABEL[v]} | {fmt(d['calibration'])} | {fmt(d['power'])} | {fmt(d['biology'])} | {d['verdict'] or 'pending'} |")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "README.md").write_text("\n".join(L) + "\n")
    print("wrote", OUT / "README.md")
