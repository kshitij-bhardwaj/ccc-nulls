"""AMENDMENT_9's report: benchmarks/results/liana_final/amendment9/README.md, from its result files and AMENDMENT_8's.

    python -m benchmarks.liana_final.amendment9 --part report
"""
from __future__ import annotations

import json

import numpy as np

from . import pipeline as P

R8, R9 = P.OUT / "amendment8", P.OUT / "amendment9"
CUTS = ("0.005", "0.01", "0.02", "0.05", "0.1")
CUT_LABEL = {"0.005": "0.5%", "0.01": "1%", "0.02": "2%", "0.05": "5%", "0.1": "10%"}
NULLS = ("null_b", "b2", "t2")
NAME = {"null_b": "BPT", "b2": "MDT", "t2": "LOO-t"}
STEPS = ("swap", "swap3", "swap4")
STEP_LABEL = {"swap": "2 tools", "swap3": "3 tools", "swap4": "4 tools"}
BASE = "liana6"
TISSUES = ("tabula_sapiens_mammary", "tabula_sapiens_liver")


def J(path):
    return json.loads(path.read_text()) if path.exists() else {"runs": {}}


def runs(f: str, dataset=None) -> dict:
    return {t: r for t, r in J(R9 / f)["runs"].items() if dataset is None or r.get("dataset") == dataset}


def scenarios(rs: dict) -> dict:
    """(dataset, side) -> its replicates."""
    out = {}
    for r in rs.values():
        out.setdefault((r["dataset"], r["scenario"]), []).append(r)
    return out


def mean_at(reps, v, c):
    return float(np.mean([r["false_share_at"][v][c] for r in reps]))


def criterion(sc: dict, v: str, c: str) -> bool:
    """AMENDMENT_1's criterion 1 at cutoff c: mean false share below liana6's in every scenario of sc."""
    return bool(sc) and all(mean_at(reps, v, c) < mean_at(reps, BASE, c) for reps in sc.values())


def verdict(sc: dict, v: str) -> str:
    ok = [criterion(sc, v, c) for c in CUTS]
    if all(ok):
        return "robust to the cutoff"
    if not ok[0]:
        return f"fails at {CUT_LABEL[CUTS[0]]}"
    first = ok.index(False)
    return f"holds up to {CUT_LABEL[CUTS[first - 1]]}; fails at {CUT_LABEL[CUTS[first]]}"


def checks() -> list[str]:
    """1% against AMENDMENT_8, tag by tag: false share, power in the top 1% and AUROC."""
    bad, n = [], 0
    for f in ("calibration.json", "confirm_calibration.json", "mouse_calibration.json"):
        old = J(R8 / f)["runs"]
        for t, r in J(R9 / f)["runs"].items():
            for v, x in r["false_share_at"].items():
                o = old.get(t, {}).get("false_share_top", {}).get(v)
                if o is None:
                    continue
                n += 1
                if abs(o - x["0.01"]) > 1e-12:
                    bad.append(f"{f} {t} {v}: {x['0.01']:.6f} here, {o:.6f} in AMENDMENT_8")
    old = J(R8 / "hard_power.json")["runs"]
    for t, r in J(R9 / "hard_power.json")["runs"].items():
        for v, d in r.items():
            if not isinstance(d, dict) or "precision_at_k" not in d or v not in old.get(t, {}):
                continue
            for c, m in d.items():
                if c == "precision_at_k":
                    continue
                o = old[t][v][c]
                n += 1
                if abs(o["power_top"] - m["power_at"]["0.01"]) > 1e-12 or not np.isclose(o["auroc"], m["auroc"], rtol=0, atol=1e-12):
                    bad.append(f"hard_power.json {t} {v} {c}: power {m['power_at']['0.01']:.3f} / AUROC {m['auroc']:.4f} here, "
                               f"{o['power_top']:.3f} / {o['auroc']:.4f} in AMENDMENT_8")
    return [f"{n:,} values compared with AMENDMENT_8; {len(bad)} differ."] + bad


def share_table(sc: dict, variants) -> list[str]:
    L = ["| Variant | " + " | ".join(CUT_LABEL[c] for c in CUTS) + " |", "|---|" + "---:|" * len(CUTS)]
    for v in variants:
        cells = []
        for c in CUTS:
            ms = [mean_at(reps, v, c) for reps in sc.values()]
            mark = "" if v == BASE else (" ✓" if criterion(sc, v, c) else " ✗")
            cells.append(f"{100 * min(ms):.0f}–{100 * max(ms):.0f}%{mark}")
        L.append(f"| {v} | " + " | ".join(cells) + " |")
    return L


def changes(sc: dict, variants, title: str) -> list[str]:
    """Where criterion 1 at another cutoff differs from the 1% result, with the scenarios that decide it."""
    L = []
    for v in variants[1:]:
        at1 = criterion(sc, v, "0.01")
        for c in CUTS:
            if c == "0.01" or criterion(sc, v, c) == at1:
                continue
            bad = [f"{d.replace('tabula_sapiens_', '')} {s} ({100 * mean_at(r, v, c):.1f}% vs liana6 {100 * mean_at(r, BASE, c):.1f}%)"
                   for (d, s), r in sc.items() if not mean_at(r, v, c) < mean_at(r, BASE, c)]
            how = "fails" if at1 else "passes"
            L.append(f"- {title}, {v} at {CUT_LABEL[c]}: {how} (at 1% it {'passes' if at1 else 'fails'})"
                     + (f"; not below liana6 in {', '.join(bad)}" if bad else "") + ".")
    return L


def write():
    cal = scenarios(runs("calibration.json"))
    lung = scenarios(runs("confirm_calibration.json"))
    mouse = scenarios(runs("mouse_calibration.json"))
    variants = [BASE] + [f"{s}_{n}" for n in NULLS for s in STEPS]
    L = ["# AMENDMENT_9 results: does the top-1% cutoff decide the result?", "",
         "Plan: [`AMENDMENT_9.md`](../../../liana_final/AMENDMENT_9.md). AMENDMENT_8's runs (same seeds, shuffles and planted "
         "pairs), scored at the top 0.5%, 1%, 2%, 5% and 10% of each ranking and by `auroc_false`, which uses no cutoff.", "",
         "## Checks against AMENDMENT_8 (the 1% column)", ""]
    L += [f"- {x}" for x in checks()[:21]]
    L += ["", "## False share of the top c: range over scenarios", "",
          "Each cell: lowest–highest mean false share (3 replicates) over the scenarios. ✓: below liana6's in every "
          "scenario (criterion 1 at that cutoff); ✗: not.", ""]
    for title, sc in (("Mammary and Liver, both sides (4 scenarios)", cal), ("Held-out Lung, both sides", lung),
                      ("Mouse brain, both sides (reported, not decided)", mouse)):
        if sc:
            L += [f"### {title}", ""] + share_table(sc, variants) + [""]
    L += ["## Verdicts for BPT", "", "| Variant | Mammary and Liver | Held-out Lung |", "|---|---|---|"]
    for s in STEPS:
        v = f"{s}_null_b"
        L.append(f"| {v} ({STEP_LABEL[s]}) | {verdict(cal, v) if cal else '–'} | {verdict(lung, v) if lung else '–'} |")
    L += ["", "For MDT and LOO-t, the ✓/✗ marks above apply the same criterion (reported, not decided).", "",
          "## What changes away from 1%", "",
          "Every case where criterion 1 at another cutoff disagrees with its 1% result. The share of false rows among "
          "all rows (which every variant reaches at 100%) is " +
          ", ".join(f"{d.replace('tabula_sapiens_', '')} {s} {100 * np.mean([r['false_share_all'] for r in reps]):.0f}%"
                    for (d, s), reps in {**cal, **lung}.items()) + ".", ""]
    ch = changes(cal, variants, "Mammary and Liver") + changes(lung, variants, "Lung")
    L += ch or ["- None."]
    L += ["",
          "## No cutoff: `auroc_false`", "",
          "P(a row that is not false outranks a false row); higher is better. Mean over replicates, per scenario.", ""]
    allsc = {**cal, **lung, **mouse}
    if allsc:
        cols = list(allsc)
        L += ["| Variant | " + " | ".join(f"{d.replace('tabula_sapiens_', '').replace('_', ' ')} {s}" for d, s in cols) +
              " | above liana6 everywhere |", "|---|" + "---:|" * (len(cols) + 1)]
        mean_auc = lambda reps, v: float(np.mean([r["auroc_false"][v] for r in reps]))
        for v in variants:
            above = "" if v == BASE else ("yes" if all(mean_auc(allsc[k], v) > mean_auc(allsc[k], BASE) for k in cols) else "no")
            L.append(f"| {v} | " + " | ".join(f"{mean_auc(allsc[k], v):.3f}" for k in cols) + f" | {above} |")
    hp = runs("hard_power.json")
    if hp:
        L += ["", "## Planted pairs (A2.2's design, f = 1.25, Mammary, 5 replicates): power in the top c, single pairs", "",
              "Reported, not decided. AUROC uses no cutoff.", "",
              "| Variant | " + " | ".join(CUT_LABEL[c] for c in CUTS) + " | AUROC |", "|---|" + "---:|" * (len(CUTS) + 1)]
        for v in variants:
            rs = [r[v]["single"] for r in hp.values() if v in r]
            if rs:
                L.append(f"| {v} | " + " | ".join(f"{100 * np.mean([x['power_at'][c] for x in rs]):.0f}%" for c in CUTS) +
                         f" | {np.mean([x['auroc'] for x in rs]):.3f} |")
    (R9 / "README.md").write_text("\n".join(L) + "\n")
    print("wrote", R9 / "README.md")
