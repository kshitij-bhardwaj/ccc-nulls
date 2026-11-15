"""benchmarks/results/liana_final/amendment3/README.md from the AMENDMENT_3 result files."""
from __future__ import annotations

import json

import numpy as np

from benchmarks.liana_nulls.stages import mean_ci

from .amendment3 import NULL_B, OUT3, SHOWN, T2, TISSUES, decision
from .pipeline import OUT

NULLS = ["null_a", "null_b", "b2", "t2"]
SHORT = {"tabula_sapiens_mammary": "Mammary", "tabula_sapiens_liver": "Liver", "tnbc_wu2021": "TNBC",
         "allen_brain_atlas": "Allen brain", "moffitt2018_hypothalamus": "hypothalamus"}


def J(path):
    return json.loads(path.read_text()) if path.exists() else None


def f3(x):
    return "–" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}"


def mci(xs):
    m, (lo, hi) = mean_ci(xs)
    return f"{m:.3f} [{lo:.3f}, {hi:.3f}]"


def table(header, rows):
    return ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)] + ["| " + " | ".join(map(str, r)) + " |" for r in rows] + [""]


def cal_rows(cal, t, s):
    return [r for r in cal["runs"].values() if r["dataset"] == t and r["scenario"] == s]


def write():
    cal, pw, hp, real, bio = (J(OUT3 / f"{n}.json") for n in ("calibration", "power", "hard_power", "real", "biology"))
    L = ["# AMENDMENT_3 results: T2 inside LIANA+'s consensus", "",
         "T2 (PR #11): each row's score against the same interaction's scores at its other cell-type pairs "
         "(leave-one-out t). Plan: [`AMENDMENT_3.md`](../../../liana_final/AMENDMENT_3.md). Every number below comes from "
         "runs that also produced liana6 and the Null B variants, so the columns are paired. Intervals: 95% t-intervals "
         "over replicates. Where each run was computed: the `machine` field of each run in the JSON files.", ""]

    L += ["## Decision (AMENDMENT_1's rule)", ""]
    d = decision()
    L += table(["Variant", "1 calibration", "2 power", "3 biology", "Verdict"],
               [[v, d[v]["calibration"], d[v]["power"], d[v]["biology"], d[v]["verdict"] or "pending"] for v in T2])
    L += ["For comparison, both Null B variants: *valid, not biologically distinguishable* (PLAN.md's results).", ""]

    if cal:
        L += ["## Calibration: false share in the top 1% (lower is better)", ""]
        rows = []
        for t in TISSUES:
            for s in ("receptor", "ligand"):
                rs = cal_rows(cal, t, s)
                if rs:
                    rows.append([f"{SHORT[t]}, {s}s shuffled", len(rs)] + [mci([r["false_share_top"][v] for r in rs]) for v in SHOWN])
        L += table(["Scenario", "reps", *SHOWN], rows)
        L += ["## p-values of the false rows (share with p ≤ 0.05; a valid test gives ≤ 0.05)", ""]
        rows = []
        for tool in ("CellPhoneDB", "CellChat"):
            rs = [r["p_calibration"][tool] for r in cal["runs"].values() if "p_calibration" in r]
            if rs:
                rows.append([tool, len(rs)] + [mci([r[n]["at_most"]["0.05"] for r in rs]) for n in NULLS] +
                            [f"{np.mean([r[n]['bh_discoveries'] for r in rs]):,.0f}" for n in NULLS])
        L += table(["Tool", "runs", *[f"p ≤ 0.05, {n}" for n in NULLS], *[f"BH 5% calls, {n}" for n in NULLS]], rows)

    for name, res, classes in (("Power (F2 design)", pw, ["single"]), ("Hard power (A2.2 design)", hp, ["single", "complex", "unequal"])):
        if not res:
            continue
        L += [f"## {name}: power in the top 1% / AUROC planted vs one-sided", ""]
        rows = []
        for f in sorted({r["effect"] for r in res["runs"].values()}):
            rs = [r for r in res["runs"].values() if r["effect"] == f]
            for c in classes:
                if c not in rs[0]["liana6"]:
                    continue
                rows.append([f"f = {f}", c, len(rs)] + [f"{np.mean([r[v][c]['power_top'] for r in rs]):.2f} / "
                                                        f"{np.nanmean([r[v][c]['auroc'] for r in rs]):.3f}" for v in SHOWN])
        L += table(["Effect", "pairs", "reps", *SHOWN], rows)

    L += ["## T2 − Null B on the same runs (reported, not judged)", "",
          "Calibration: lower is better, so a negative difference favours T2. Power and AUROC: a positive difference favours T2.", ""]
    rows = []
    for (t2, nb) in zip(T2, NULL_B):
        if cal:
            for t in TISSUES:
                for s in ("receptor", "ligand"):
                    rs = cal_rows(cal, t, s)
                    if rs:
                        rows.append([f"{t2} − {nb}", f"false share top 1%, {SHORT[t]} {s}",
                                     mci([r["false_share_top"][t2] - r["false_share_top"][nb] for r in rs])])
        for name, res in (("F2", pw), ("A2.2", hp)):
            if not res:
                continue
            for f in sorted({r["effect"] for r in res["runs"].values()}):
                rs = [r for r in res["runs"].values() if r["effect"] == f]
                for m in ("power_top", "auroc"):
                    rows.append([f"{t2} − {nb}", f"{m}, {name} single, f = {f}", mci([r[t2]["single"][m] - r[nb]["single"][m] for r in rs])])
    L += table(["Difference", "Measure", "Mean [95% CI]"], rows)

    if real:
        L += ["## The five data sets: calls (specificity_rank ≤ 0.05) and Jaccard with liana6", ""]
        L += table(["Data set", *SHOWN, "T2 untestable rows (CellPhoneDB / CellChat)"],
                   [[SHORT[n]] + [f"{r['calls'][v]['calls']:,} ({r['calls'][v]['jaccard_with_liana6']:.2f})" for v in SHOWN] +
                    [" / ".join(f"{r['tools'][t]['t2_untested']:.1%}" for t in ("CellPhoneDB", "CellChat"))]
                    for n, r in real["runs"].items()])
    if bio:
        L += ["## Biology (F3): spatial odds ratio, liana6 calls kept vs dropped [95% CI]", ""]
        L += table(["Data set", *SHOWN[1:]],
                   [[SHORT[n]] + [("–" if x.get("odds_ratio") is None else f"{x['odds_ratio']:.2f} [{x['ci95'][0]:.2f}, {x['ci95'][1]:.2f}]")
                                  for x in (bio["spatial"][n][v] for v in SHOWN[1:])] for n in bio["spatial"]])
        L += ["Downstream response, odds ratio kept vs dropped:", ""]
        L += table(["Data set", "footprints", *SHOWN[1:]],
                   [[SHORT[n], res] + [f3(bio["downstream"][n][res][v].get("odds_ratio")) for v in SHOWN[1:]]
                    for n in bio["downstream"] for res in ("cytosig", "progeny")])

    L += ["## Checks: the earlier variants reproduce", ""]
    rows = []
    rc, hp0, pw0 = J(OUT / "amendment2" / "recalibration.json"), J(OUT / "amendment2" / "hard_power.json"), J(OUT / "power.json")
    for label, new, old, get in (
            ("calibration vs A2.6", cal, rc, lambda r, v: r["false_share_top"][v]),
            ("hard power vs A2.2", hp, hp0, lambda r, v: (r[v]["single"]["power_top"], r[v]["single"]["auroc"]))):
        if new and old:
            both = [t for t in new["runs"] if t in old["runs"]]
            same = sum(all(get(new["runs"][t], v) == get(old["runs"][t], v) for v in ["liana6", *NULL_B]) for t in both)
            rows.append([label, f"{same} of {len(both)} runs identical"])
    if pw and pw0:
        diffs = []
        for f in ("1.5", "2.0"):
            rs = [r for r in pw["runs"].values() if str(r["effect"]) == f]
            for v in ["liana6", *NULL_B]:
                if rs and f in pw0["summary"]:
                    diffs.append(abs(np.mean([r[v]["single"]["power_top"] for r in rs]) - pw0["summary"][f][v]["power_top"][0]))
        rows.append(["power vs F2 (means)", f"largest difference {max(diffs):.4f}" if diffs else "–"])
    if real:
        rows += [[f"real {SHORT[n]} vs run tables", r.get("same_as_run", "run table not available")] for n, r in real["runs"].items()]
    L += table(["Check", "Result"], rows)
    OUT3.mkdir(parents=True, exist_ok=True)
    (OUT3 / "README.md").write_text("\n".join(L) + "\n")
    print("wrote", OUT3 / "README.md")
