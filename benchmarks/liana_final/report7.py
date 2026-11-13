"""benchmarks/results/liana_final/amendment7/README.md from the AMENDMENT_7 result files.

    python -m benchmarks.liana_final.report7
"""
from __future__ import annotations

import json

import numpy as np

from .amendment2 import bh
from .amendment6 import SHOWN
from .amendment7 import OUT7, TOOL_NULLS
from .pipeline import BASELINE, OUT
from .report3 import mci, table

OUT2, OUT4, OUT5 = OUT / "amendment2", OUT / "amendment4", OUT / "amendment5"
TOOLS = ("CellPhoneDB", "CellChat", "Connectome")
SHORT = {"tabula_sapiens_lung": "Lung", "tnbc_wu2021": "TNBC", "tabula_sapiens_mammary": "Mammary",
         "tabula_sapiens_liver": "Liver", "allen_brain_atlas": "Allen brain", "moffitt2018_hypothalamus": "hypothalamus"}
NAME = {"null_a": "JLP (#/B)", "null_a_plus1": "JLP ((1+#)/(B+1))", "null_b": "BPT", "b2": "MDT", "t2": "LOO-t"}


def J(path):
    return json.loads(path.read_text()) if path.exists() else None


def f3(x):
    return "–" if x is None else f"{x:.3f}"


def donor(L):
    dn = J(OUT7 / "donor_null.json")
    if not dn:
        return
    L += ["## Donors (A2.3 a): labels shuffled within donor; share of rows with p ≤ 0.05", "",
          "s ≠ t, magnitude > 0. Global: labels permuted over all cells (LIANA+'s scheme); donor: within each donor. "
          "MDT uses no permutations, so it has the global column only.", ""]
    rows = []
    for name in dict.fromkeys(r["dataset"] for r in dn["runs"].values()):
        for tool in TOOLS:
            g = [r for r in dn["runs"].values() if r["dataset"] == name and r["scheme"] == "global"]
            d = [r for r in dn["runs"].values() if r["dataset"] == name and r["scheme"] == "donor"]
            rows.append([SHORT[name], tool, f"{len(g)} / {len(d)}"] +
                        [f"{mci([r[tool][f'{n}_le_0.05'] for r in g])} / {mci([r[tool][f'{n}_le_0.05'] for r in d])}"
                         for n in ("null_a", "null_b", "t2")] + [mci([r[tool]["b2_le_0.05"] for r in g])])
    L += table(["Data set", "Tool", "reps (global / donor)", "JLP: global / donor", "BPT: global / donor",
                "LOO-t: global / donor", "MDT"], rows)
    return dn


def donor_checks(dn):
    """CellPhoneDB's and CellChat's numbers against AMENDMENT_2 (JLP, BPT), AMENDMENT_4 (LOO-t), AMENDMENT_5 (MDT)."""
    out = []
    for label, path, keys, tag in (
            ("AMENDMENT_2 (JLP, BPT; ≤ 0.01 and ≤ 0.05)", OUT2 / "donor_null.json",
             [f"{n}_le_{a}" for n in ("null_a", "null_b") for a in (0.01, 0.05)], lambda t: t),
            ("AMENDMENT_4 (LOO-t)", OUT4 / "donor_null.json", ["t2_le_0.01", "t2_le_0.05"], lambda t: t),
            ("AMENDMENT_5 (MDT)", OUT5 / "donor_b2.json", ["b2_le_0.01", "b2_le_0.05"],
             lambda t: t.replace("|global|", "|") if "|global|" in t else None)):
        old = J(path)
        if not old:
            out.append([label, "file missing"])
            continue
        pairs = [(t, tag(t)) for t in dn["runs"] if tag(t) in old["runs"]]
        same = sum(all(dn["runs"][t][tool].get(k) == old["runs"][o][tool].get(k) for tool in ("CellPhoneDB", "CellChat")
                       for k in keys if k in old["runs"][o][tool]) for t, o in pairs)
        out.append([label, f"{same} of {len(pairs)} runs identical"])
    return out


def sensitivity(L):
    f = lambda n: J(OUT7 / f"{n}.json")
    sc, sp, ss, sP, sPc = (f(n) for n in ("sensitivity_calibration", "sensitivity_power", "sensitivity_seed",
                                          "sensitivity_permutations", "sensitivity_permutations_calibration"))
    if not (sc or sp or ss or sP or sPc):
        return
    L += ["## Sensitivity (A2.5), with the swap3 variants", "",
          "✓: false share of the top 1% below liana6's in the same run; ✗: not below.", ""]
    for name, res in (("Calibration: false share in the top 1%", sc), ("Permutations, calibration (receptors, rep 0)", sPc)):
        if res:
            L += [f"{name}:", ""]
            L += table(["Run", *SHOWN], [[t] + [f"{r['false_share_top'][v]:.3f}" + ("" if v == BASELINE else
                                                 " ✓" if r["false_share_top"][v] < r["false_share_top"][BASELINE] else " ✗")
                                                 for v in SHOWN] for t, r in res["runs"].items()])
    for name, res in (("Calibration", sc), ("Permutations, calibration", sPc)):
        if res:
            n = len(res["runs"])
            L += [f"{name}: runs where the variant's false share is below liana6's — " + ", ".join(
                f"{v} {sum(r['false_share_top'][v] < r['false_share_top'][BASELINE] for r in res['runs'].values())} of {n}"
                for v in SHOWN if v != BASELINE) + ".", ""]
    if sp:
        L += ["Power, f = 1.5 (power in the top 1% / AUROC, single pairs):", ""]
        L += table(["Run", *SHOWN], [[t] + [f"{r[v]['single']['power_top']:.2f} / {r[v]['single']['auroc']:.3f}" for v in SHOWN]
                                    for t, r in sp["runs"].items()])
    for name, res in (("Seeds, real labels (hypothalamus), against seed 1337", ss),
                      ("5,000 vs 1,000 permutations, real labels", sP)):
        if res:
            L += [f"{name} (Spearman of specificity_rank / top-1% Jaccard):", ""]
            L += table(["Run", *SHOWN], [[t] + [f"{r[v]['spearman']:.4f} / {r[v]['top1_jaccard']:.3f}" if v in r else "not run"
                                                for v in SHOWN] for t, r in res["runs"].items()])
    if ss or sP:
        L += ["\"not run\": A2.5's seed and permutation runs compute no MDT (`b2=False`, as in AMENDMENT_4), so `swap_b2` "
              "and `swap3_b2` have no second table to compare. The first seed run stopped on this (`KeyError`); "
              "`stability` was changed to record these variants as skipped (commit 42b952f), and the run was repeated. "
              "No other result changed: the sessions that had finished ran the same code for every part they computed.", ""]

    def same(new, old, get):
        """Runs in both files whose earlier variants (those the old run has) are identical."""
        if not (new and old):
            return None
        both = [t for t in new["runs"] if t in old["runs"]]
        n = sum(all(get(new["runs"][t], v) == get(old["runs"][t], v) for v in early if get(old["runs"][t], v) is not None)
                for t in both)
        return f"{n} of {len(both)} runs identical"
    early = [BASELINE, "swap_null_b", "swap_b2", "swap_t2"]
    fs = lambda r, v: r["false_share_top"].get(v)
    pw = lambda r, v: (r[v]["single"]["power_top"], r[v]["single"]["auroc"]) if v in r else None
    st = lambda r, v: (r[v]["spearman"], r[v]["top1_jaccard"]) if v in r else None
    rows = [[n, same(res, J(OUT4 / f"{n}.json"), g) or "–"] for n, res, g in (
        ("sensitivity_calibration", sc, fs), ("sensitivity_power", sp, pw), ("sensitivity_seed", ss, st),
        ("sensitivity_permutations", sP, st), ("sensitivity_permutations_calibration", sPc, fs))]
    return rows


def cases(L):
    cs = J(OUT7 / "cases.json")
    if not cs:
        return
    L += ["## Case studies (A2.9): what swap3_<n> drops among liana6's calls (s ≠ t)", "",
          "Breadth: share of cell types where the gene is detected in ≥ 10% of cells; a row's breadth is the larger of its "
          "ligand's and receptor's. p: one-sided Mann–Whitney, dropped broader than kept. Overlap with what swap_<n> drops "
          "in the same table. Computed locally from AMENDMENT_6's saved real tables, not in a Kaggle session "
          "(`python -m benchmarks.liana_final.amendment7 --part cases`).", ""]
    rows = []
    for name, r in cs["runs"].items():
        for v in ("swap3_null_b", "swap3_b2", "swap3_t2"):
            x = r[v]
            rows.append([SHORT[name], v, x["kept"], x["dropped"],
                         f"{f3(x['breadth_dropped_median'])} / {f3(x['breadth_kept_median'])}",
                         "–" if x["mann_whitney_p"] is None else f"{x['mann_whitney_p']:.2g}",
                         f"{x['dropped_by_both']} / {x['dropped_swap3_only']} / {x['dropped_swap_only']}"])
    L += table(["Data set", "Variant", "kept", "dropped", "breadth dropped / kept", "p",
                "dropped by both / swap3 only / swap only"], rows)


def tool_level(L):
    tl, tp, tr = (J(OUT7 / f"{n}.json") for n in ("tool_labels", "tool_power", "tool_real"))
    if tl or tp or tr:
        L += ["## Benchmark 1's tests at tool level, Connectome beside CellPhoneDB and CellChat", "",
              "Each tool's own rows inside LIANA+, p from each null, BH at 5% over the rows with a p. JLP twice: LIANA+'s "
              "p (#/B, can be 0) and Benchmark 1's ((1 + #)/(B + 1)).", ""]
    if tl:
        L += ["### Labels shuffled (every label permuted; any call is a false call)", "",
              "Runs with any call / mean number of calls.", ""]
        rows = []
        for name in dict.fromkeys(r["dataset"] for r in tl["runs"].values()):
            rs = [r for r in tl["runs"].values() if r["dataset"] == name]
            for tool in TOOLS:
                rows.append([SHORT[name], tool, len(rs)] + [f"{sum(r[tool][n]['any_call'] for r in rs)} / "
                                                            f"{np.mean([r[tool][n]['calls'] for r in rs]):,.1f}" for n in TOOL_NULLS])
        L += table(["Data set", "Tool", "reps", *(NAME[n] for n in TOOL_NULLS)], rows)
    if tp:
        L += ["### Power (F2's design, 30 single pairs, f = 1.5, Mammary): power / FDP, mean over replicates", ""]
        rows = []
        for B in sorted({r["n_perms"] for r in tp["runs"].values()}):
            rs = [r for r in tp["runs"].values() if r["n_perms"] == B]
            for tool in TOOLS:
                rows.append([f"{B:,}", tool, len(rs)] + [f"{np.mean([r[tool][n]['power'] for r in rs]):.2f} / "
                                                         f"{np.mean([r[tool][n]['fdp'] for r in rs]):.2f}" for n in TOOL_NULLS])
        L += table(["B", "Tool", "reps", *(NAME[n] for n in TOOL_NULLS)], rows)
    if tr:
        L += ["### Real labels: BH calls", ""]
        L += table(["Data set", "Tool", "rows", *(NAME[n] for n in TOOL_NULLS)],
                   [[SHORT[name], tool, f"{r[tool]['rows']:,}"] + [f"{r[tool][f'calls_{n}']:,}" for n in TOOL_NULLS]
                    for name, r in tr["runs"].items() for tool in TOOLS if tool in r])
        # each session wrote q over its own data set only; recomputed here over every spatial test of every session
        tests = [(k, t, n) for k, r in tr["runs"].items() for t in TOOLS if "spatial" in r.get(t, {})
                 for n in r[t]["spatial"] if r[t]["spatial"][n].get("odds_ratio") is not None]
        q = dict(zip(tests, bh(np.array([tr["runs"][k][t]["spatial"][n]["p"] for k, t, n in tests])))) if tests else {}
        if tests:
            L += ["### Spatial truth: among JLP's calls ((1 + #)/(B + 1)), odds ratio of the rows each null keeps vs drops "
                  "[95% CI]; q: BH over all tests below", ""]
            rows = []
            for k, r in tr["runs"].items():
                for t in TOOLS:
                    if "spatial" not in r.get(t, {}):
                        continue
                    rows.append([SHORT[k], t] + [(lambda x, key: f"{x['odds_ratio']:.2f} [{x['ci95'][0]:.2f}, {x['ci95'][1]:.2f}] "
                                                  f"q {q[key]:.2f}" if x.get("odds_ratio") is not None else "–")(
                        r[t]["spatial"].get(n, {}), (k, t, n)) for n in TOOL_NULLS[2:]])
            L += table(["Data set", "Tool", *(NAME[n] for n in TOOL_NULLS[2:])], rows)
            sig = [f"{SHORT[k]} {t} {NAME[n]}" for (k, t, n), x in q.items() if x < 0.05]
            L += [f"Tests with q < 0.05: {', '.join(sig) if sig else 'none'}.", ""]


def write():
    L = ["# AMENDMENT_7 results: Connectome through the remaining tests", "",
         "Plan: [`AMENDMENT_7.md`](../../../liana_final/AMENDMENT_7.md). Reported, not decided: AMENDMENT_6's decision for "
         "`swap3_<n>` stands. Names as in the paper: JLP = Null A (LIANA+'s own test), BPT = Null B, MDT = B2, "
         "LOO-t = T2. Intervals: 95% t-intervals over replicates.", ""]
    dn = donor(L)
    sens_checks = sensitivity(L)
    cases(L)
    tool_level(L)
    L += ["## Checks: CellPhoneDB's and CellChat's numbers, and the earlier variants, reproduce earlier amendments", ""]
    rows = (donor_checks(dn) if dn else []) + [[f"AMENDMENT_4 {n} (liana6, swap_null_b, swap_b2, swap_t2)", r]
                                               for n, r in (sens_checks or [])]
    L += table(["Against", "Result"], rows)
    OUT7.mkdir(parents=True, exist_ok=True)
    (OUT7 / "README.md").write_text("\n".join(L) + "\n")
    print("wrote", OUT7 / "README.md")


if __name__ == "__main__":
    write()
