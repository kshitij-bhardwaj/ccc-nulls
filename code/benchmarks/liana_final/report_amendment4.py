"""benchmarks/results/liana_final/amendment4/README.md from the AMENDMENT_4 result files."""
from __future__ import annotations

import json

import numpy as np

from .amendment4 import NULL_B, OUT4, SHOWN, T2
from .pipeline import BASELINE, OUT
from .report_amendment3 import SHORT, f3, mci, table

SHORT = {**SHORT, "tabula_sapiens_lung": "Lung"}
OUT2 = OUT / "amendment2"


def J(path):
    return json.loads(path.read_text()) if path.exists() else None


def crit1(cal, name, v, setting=None):
    """AMENDMENT_1's criterion 1 on one data set: below liana6 on both sides (None while runs are missing)."""
    ok = []
    for s in ("receptor", "ligand"):
        rs = [r for r in cal["runs"].values() if r["dataset"] == name and r["scenario"] == s
              and (setting is None or r.get("setting") == setting)]
        if len(rs) < 3:
            return None
        ok.append(np.mean([r["false_share_top"][v] for r in rs]) < np.mean([r["false_share_top"][BASELINE] for r in rs]))
    return all(ok)


def crit2(pw, v):
    rs = [r for r in pw["runs"].values() if r["effect"] == 2.0]
    if len(rs) < 5:
        return None
    return bool(np.mean([r[v]["single"]["power_top"] for r in rs]) >= 0.8 and
                np.mean([r[v]["single"]["auroc"] for r in rs]) >= np.mean([r[BASELINE]["single"]["auroc"] for r in rs]))


def cal_table(cal, name):
    rows = []
    for s in ("receptor", "ligand"):
        rs = [r for r in cal["runs"].values() if r["dataset"] == name and r["scenario"] == s and r.get("setting", "default") == "default"]
        if rs:
            rows.append([f"{s}s shuffled", len(rs)] + [mci([r["false_share_top"][v] for r in rs]) for v in SHOWN])
    return table(["Scenario", "reps", *SHOWN], rows)


def pw_table(pw):
    rows = []
    for f in sorted({r["effect"] for r in pw["runs"].values()}):
        rs = [r for r in pw["runs"].values() if r["effect"] == f]
        rows.append([f"f = {f}", len(rs)] + [f"{np.mean([r[v]['single']['power_top'] for r in rs]):.2f} / "
                                             f"{np.nanmean([r[v]['single']['auroc'] for r in rs]):.3f}" for v in SHOWN])
    return table(["Effect", "reps", *SHOWN], rows)


def word(x, yes, no):
    return "pending" if x is None else yes if x else no


def same(new, old, get):
    """Runs in both files whose liana6 and Null B values are identical."""
    if not (new and old):
        return None
    both = [t for t in new["runs"] if t in old["runs"]]
    try:
        n = sum(all(get(new["runs"][t], v) == get(old["runs"][t], v) for v in [BASELINE, *NULL_B]) for t in both)
    except (KeyError, TypeError):
        return "not comparable"
    return f"{n} of {len(both)} runs identical"


def write():
    f = lambda n: J(OUT4 / f"{n}.json")
    L = ["# AMENDMENT_4 results: T2 on the rest of AMENDMENT_2", "",
         "Plan: [`AMENDMENT_4.md`](../../../liana_final/AMENDMENT_4.md). T2 as in AMENDMENT_3; every run also produced "
         "liana6 and the Null B variants, so the columns are paired. Intervals: 95% t-intervals over replicates. "
         "Power cells: power in the top 1% / AUROC planted vs one-sided.", ""]

    cc, cp = f("confirm_calibration"), f("confirm_power")
    L += ["## Lung (A2.1): held-out confirmation", ""]
    rows = []
    for v in [*NULL_B, *T2]:
        c1 = crit1(cc, "tabula_sapiens_lung", v) if cc else None
        c2 = crit2(cp, v) if cp else None
        out = None if None in (c1, c2) else "confirmed" if c1 and c2 else \
            "not confirmed (" + ", ".join(n for n, c in (("calibration", c1), ("power", c2)) if not c) + ")"
        rows.append([v, word(c1, "holds", "fails"), word(c2, "holds", "fails"), out or "pending"])
    L += table(["Variant", "1 calibration", "2 power", "Outcome"], rows)
    if cc:
        L += cal_table(cc, "tabula_sapiens_lung")
    if cp:
        L += pw_table(cp)

    dn, dr = f("donor_null"), f("donor_real")
    if dn:
        L += ["## Donors (A2.3 a): labels shuffled within donor, share of rows with p ≤ 0.05 (nominal 0.05)", ""]
        rows = []
        for name in dict.fromkeys(r["dataset"] for r in dn["runs"].values()):
            rs = [r for r in dn["runs"].values() if r["dataset"] == name]
            for tool in ("CellPhoneDB", "CellChat"):
                rows.append([SHORT.get(name, name), tool, len(rs)] +
                            [mci([r[tool][f"{n}_le_0.05"] for r in rs]) for n in ("null_a", "null_b", "t2")] +
                            [f"{np.mean([r[tool]['t2_untested'] for r in rs]):.1%}"])
        L += table(["Data set", "Tool", "reps", "Null A", "Null B", "T2", "T2 untestable"], rows)
    if dr:
        L += ["## Donors (A2.3 b): real labels, within-donor vs global permutations", ""]
        L += table(["Data set", *[f"{v}: top-1% / call Jaccard" for v in ["swap_null_b", *T2]]],
                   [[SHORT[n]] + [f"{r[v]['top1_jaccard']:.3f} / {r[v]['call_jaccard']:.3f}" for v in ["swap_null_b", *T2]]
                    for n, r in dr["runs"].items()])

    sc, sp, ss, sP, sPc = (f(n) for n in ("sensitivity_calibration", "sensitivity_power", "sensitivity_seed",
                                          "sensitivity_permutations", "sensitivity_permutations_calibration"))
    if sc or sp or ss or sP:
        L += ["## Sensitivity (A2.5)", ""]
    for name, res in (("Calibration: false share in the top 1% (below liana6?)", sc),
                      ("Permutations, calibration (receptors, rep 0)", sPc)):
        if res:
            L += [f"{name}:", ""]
            L += table(["Run", *SHOWN], [[t] + [f"{r['false_share_top'][v]:.3f}" + ("" if v == BASELINE else
                                                 " ✓" if r["false_share_top"][v] < r["false_share_top"][BASELINE] else " ✗")
                                                 for v in SHOWN] for t, r in res["runs"].items()])
    if sp:
        L += ["Power, f = 1.5:", ""]
        L += table(["Run", *SHOWN], [[t] + [f"{r[v]['single']['power_top']:.2f} / {r[v]['single']['auroc']:.3f}" for v in SHOWN]
                                    for t, r in sp["runs"].items()])
    for name, res in (("Seeds, real labels, against seed 1337", ss), ("5,000 vs 1,000 permutations, real labels", sP)):
        if res:
            L += [f"{name} (Spearman of specificity_rank / top-1% Jaccard):", ""]
            L += table(["Run", *SHOWN], [[t] + [f"{r[v]['spearman']:.4f} / {r[v]['top1_jaccard']:.3f}" for v in SHOWN]
                                        for t, r in res["runs"].items()])

    mc = f("mouse_calibration")
    if mc:
        L += ["## Mouse (A2.7): Allen brain, criterion 1", ""]
        L += table(["Variant", "criterion 1"], [[v, word(crit1(mc, "allen_brain_atlas", v), "holds", "fails")] for v in [*NULL_B, *T2]])
        L += cal_table(mc, "allen_brain_atlas")

    cs = f("cases")
    if cs:
        L += ["## What swap_t2 drops (among liana6's calls, s ≠ t)", "",
              "Breadth: share of cell types where the gene is detected in ≥ 10% of cells; a row's breadth is the larger "
              "of its ligand's and receptor's. p: one-sided Mann–Whitney, dropped broader than kept.", ""]
        L += table(["Data set", "kept", "dropped", "breadth dropped / kept", "p", "dropped by both / T2 only / Null B only",
                    "Jaccard of dropped sets", "promoted"],
                   [[SHORT[n], r["kept"], r["dropped"], f"{f3(r['breadth_dropped_median'])} / {f3(r['breadth_kept_median'])}",
                     "–" if r["mann_whitney_p"] is None else f"{r['mann_whitney_p']:.2g}",
                     f"{r['dropped_by_both']} / {r['dropped_t2_only']} / {r['dropped_null_b_only']}", f3(r["dropped_jaccard"]),
                     r["promoted_rows"]] for n, r in cs["runs"].items()])
        L += ["Tables: `cases_<data set>_dropped.csv` and `cases_<data set>_promoted.csv` in this folder.", ""]

    tm, tm2 = f("timing"), J(OUT2 / "timing.json")
    if tm:
        L += ["## Timing", ""]
        L += table(["Data set", "T2 (s)", "this pipeline, A2.10 (s)", "rows (CellPhoneDB / CellChat)"],
                   [[SHORT.get(n, n), r["t2_seconds"], (tm2 or {}).get("runs", {}).get(n, {}).get("this_pipeline", "–"),
                     " / ".join(f"{x:,}" for x in r["rows"].values())] for n, r in tm["runs"].items()])

    L += ["## Checks: liana6 and the Null B variants reproduce AMENDMENT_2", ""]
    fs = lambda r, v: r["false_share_top"][v]
    pw = lambda r, v: (r[v]["single"]["power_top"], r[v]["single"]["auroc"])
    rows = [[n, same(f(n), J(OUT2 / f"{n}.json"), g) or "–"] for n, g in (
        ("confirm_calibration", fs), ("confirm_power", pw), ("mouse_calibration", fs), ("sensitivity_calibration", fs),
        ("sensitivity_power", pw), ("sensitivity_permutations_calibration", fs))]
    dn2 = J(OUT2 / "donor_null.json")
    if dn and dn2:
        both = [t for t in dn["runs"] if t in dn2["runs"]]
        n = sum(all(dn["runs"][t][tool][k] == dn2["runs"][t][tool][k] for tool in ("CellPhoneDB", "CellChat")
                    for k in ("null_a_le_0.05", "null_b_le_0.05", "null_a_le_0.01", "null_b_le_0.01")) for t in both)
        rows.append(["donor_null (Null A, Null B, global)", f"{n} of {len(both)} runs identical"])
    L += table(["File", "Result"], rows)
    OUT4.mkdir(parents=True, exist_ok=True)
    (OUT4 / "README.md").write_text("\n".join(L) + "\n")
    print("wrote", OUT4 / "README.md")
