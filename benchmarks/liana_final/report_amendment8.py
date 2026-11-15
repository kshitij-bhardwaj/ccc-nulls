"""benchmarks/results/liana_final/amendment8/README.md from the AMENDMENT_8 result files.

    python -m benchmarks.liana_final.report_amendment8
"""
from __future__ import annotations

import json

import numpy as np

from .amendment2 import bh
from .amendment7 import TOOL_NULLS
from .amendment8 import NULLS4, OUT8, SWAP3, SWAP4
from .pipeline import BASELINE, OUT
from .report_amendment3 import mci, table
from .report_amendment4 import crit1, crit2, word

OUT6, OUT7 = OUT / "amendment6", OUT / "amendment7"
TOOLS = ("CellPhoneDB", "CellChat", "Connectome", "NATMI")
SWAP2 = [f"swap_{n}" for n in NULLS4 if n != "null_a"]
EARLIER = [BASELINE, *SWAP2, *SWAP3]                     # must equal AMENDMENT_6's and AMENDMENT_7's
SHOWN = [BASELINE, *[v for n in ("null_b", "b2", "t2") for v in (f"swap3_{n}", f"swap4_{n}")], "swap4_null_a"]
SHORT = {"tabula_sapiens_lung": "Lung", "tnbc_wu2021": "TNBC", "tabula_sapiens_mammary": "Mammary",
         "tabula_sapiens_liver": "Liver", "allen_brain_atlas": "Allen brain", "moffitt2018_hypothalamus": "hypothalamus"}
NAME = {"null_a": "JLP (#/B)", "null_a_plus1": "JLP ((1+#)/(B+1))", "null_b": "BPT", "b2": "MDT", "t2": "LOO-t"}
SETTINGS = [("calibration", "tabula_sapiens_mammary", "Mammary"), ("calibration", "tabula_sapiens_liver", "Liver"),
            ("confirm_calibration", "tabula_sapiens_lung", "Lung (held out)"), ("mouse_calibration", "allen_brain_atlas", "Allen brain")]


def J(path):
    return json.loads(path.read_text()) if path.exists() else None


f8 = lambda n: J(OUT8 / f"{n}.json")


def same_runs(new, old, get, variants):
    """Runs in both files whose given variants are identical (only the variants the old run has)."""
    if not (new and old):
        return "pending"
    both = [t for t in new["runs"] if t in old["runs"]]
    n = sum(all(get(new["runs"][t], v) == get(old["runs"][t], v) for v in variants if get(old["runs"][t], v) is not None)
            for t in both)
    return f"{n} of {len(both)} runs identical"


def consensus(L):
    L += ["## Calibration: false share of the top 1% (lower is better)", "",
          "swap3: CellPhoneDB's, CellChat's and Connectome's p replaced; swap4: NATMI's too. Same runs, so paired.", ""]
    rows = []
    for part, name, short in SETTINGS:
        cal = f8(part)
        for s in ("receptor", "ligand"):
            rs = [r for r in (cal or {"runs": {}})["runs"].values() if r["dataset"] == name and r["scenario"] == s
                  and r.get("setting", "default") == "default"]
            if rs:
                rows.append([f"{short}, {s}s", len(rs)] + [mci([r["false_share_top"][v] for r in rs]) for v in SHOWN])
    L += table(["Setting", "reps", *SHOWN], rows)

    cal, pw, cc, cp = f8("calibration"), f8("power"), f8("confirm_calibration"), f8("confirm_power")
    L += ["## AMENDMENT_1's rule for swap4 (and swap3, from the same runs)", "",
          "Criterion 1: false share of the top 1% below liana6's on both sides; criterion 2: power in the top 1% ≥ 0.8 "
          "and AUROC ≥ liana6's at f = 2.0.", ""]
    rows = []
    for v in [*SWAP4, *SWAP3]:
        c1 = [crit1(cal, d, v) if cal else None for d in ("tabula_sapiens_mammary", "tabula_sapiens_liver")]
        c2 = crit2(pw, v) if pw else None
        l1 = crit1(cc, "tabula_sapiens_lung", v) if cc else None
        l2 = crit2(cp, v) if cp else None
        out = None if None in (l1, l2) else "confirmed" if l1 and l2 else \
            "not confirmed (" + ", ".join(n for n, c in (("calibration", l1), ("power", l2)) if not c) + ")"
        rows.append([v, " / ".join(word(c, "holds", "fails") for c in c1), word(c2, "holds", "fails"),
                     word(l1, "holds", "fails"), word(l2, "holds", "fails"), out or "pending"])
    L += table(["Variant", "1, Mammary / Liver", "2, Mammary", "1, Lung", "2, Lung", "Held-out Lung"], rows)

    for part, title in (("power", "Power, F2 design"), ("hard_power", "Power, small effects (A2.2)"), ("confirm_power", "Power, held-out Lung")):
        p = f8(part)
        if not p:
            continue
        L += [f"## {title}: power in the top 1% / AUROC, single pairs", ""]
        rows = []
        for e in sorted({r["effect"] for r in p["runs"].values()}):
            rs = [r for r in p["runs"].values() if r["effect"] == e]
            rows.append([f"f = {e}", len(rs)] + [f"{np.mean([r[v]['single']['power_top'] for r in rs]):.2f} / "
                                                f"{np.nanmean([r[v]['single']['auroc'] for r in rs]):.3f}" for v in SHOWN])
        L += table(["Effect", "reps", *SHOWN], rows)

    if cal:
        L += ["## Each tool's p inside LIANA+ on false rows (Mammary and Liver, both sides): share ≤ 0.05 / BH discoveries / "
              "false share of discoveries", ""]
        rs = [r for r in cal["runs"].values() if "p_calibration" in r]
        rows = []
        for tool in ("Connectome", "NATMI"):
            for n in NULLS4:
                x = [r["p_calibration"][tool][n] for r in rs if tool in r["p_calibration"]]
                fs = [y["bh_false_share"] for y in x if y["bh_false_share"] is not None]
                rows.append([tool, NAME[n], len(x), mci([y["at_most"]["0.05"] for y in x]),
                             f"{np.mean([y['bh_discoveries'] for y in x]):,.0f}", f"{np.mean(fs):.3f}" if fs else "–"])
        L += table(["Tool", "Null", "runs", "p ≤ 0.05", "BH discoveries", "false share of discoveries"], rows)

    real = f8("real")
    if real:
        L += ["## Real labels: calls (specificity_rank ≤ 0.05) and Jaccard with liana6", ""]
        L += table(["Data set", *SHOWN], [[SHORT[n]] + [f"{r['calls'][v]['calls']:,} ({r['calls'][v]['jaccard_with_liana6']:.2f})"
                                                        for v in SHOWN] for n, r in real["runs"].items()])
    bio = f8("biology")
    if bio:
        t = bio["tests"]
        L += ["## Biology (A2.8): spatial odds ratio, kept vs dropped [95% CI], q across all", ""]
        rows = []
        for v in ["swap_null_b", "swap_b2", "swap_t2", *SWAP3, *SWAP4]:
            rows.append([v] + [(lambda r: f"{r['odds_ratio']:.2f} [{r['ci95'][0]:.2f}, {r['ci95'][1]:.2f}] q {r['q_all']:.2f}"
                                if r.get("odds_ratio") is not None else "–")(t.get(f"{v} | {d} spatial", {}))
                               for d in ("tnbc_wu2021", "allen_brain_atlas", "moffitt2018_hypothalamus")])
        L += table(["Variant", "TNBC", "Allen brain", "hypothalamus"], rows)
        sig = [k for k, r in t.items() if r.get("q_all") is not None and r["q_all"] < 0.05]
        L += [f"Tests with q (all) < 0.05: {', '.join(sig) if sig else 'none'}.", ""]
        down = [(v, k.split(" | ")[1], r) for k, r in t.items() for v in [k.split(" | ")[0]] if v in SWAP4
                and "spatial" not in k and r.get("odds_ratio") is not None]
        if down:
            fav = [f"{v} {d} ({r['odds_ratio']:.2f}, p {r['p']:.2g})" for v, d, r in down if r["p"] < 0.05]
            L += [f"Downstream tests (CytoSig, PROGENy) for swap4 with bootstrap p < 0.05: {', '.join(fav) if fav else 'none'}.", ""]


def donor(L):
    dn = f8("donor_null")
    if not dn:
        return
    L += ["## Donors (A2.3 a): labels shuffled within donor; rows with p ≤ 0.05", "",
          "Global: labels permuted over all cells (LIANA+'s scheme); donor: within each donor. MDT uses no permutations.", ""]
    rows = []
    for name in dict.fromkeys(r["dataset"] for r in dn["runs"].values()):
        g = [r for r in dn["runs"].values() if r["dataset"] == name and r["scheme"] == "global"]
        d = [r for r in dn["runs"].values() if r["dataset"] == name and r["scheme"] == "donor"]
        for tool in TOOLS:
            rows.append([SHORT[name], tool, f"{len(g)} / {len(d)}"] +
                        [f"{mci([r[tool][f'{n}_le_0.05'] for r in g])} / {mci([r[tool][f'{n}_le_0.05'] for r in d])}"
                         for n in ("null_a", "null_b", "t2")] + [mci([r[tool]["b2_le_0.05"] for r in g])])
    L += table(["Data set", "Tool", "reps (global / donor)", "JLP: global / donor", "BPT: global / donor",
                "LOO-t: global / donor", "MDT"], rows)


def sensitivity(L):
    sc, sp, ss, sP, sPc = (f8(n) for n in ("sensitivity_calibration", "sensitivity_power", "sensitivity_seed",
                                           "sensitivity_permutations", "sensitivity_permutations_calibration"))
    if not any((sc, sp, ss, sP, sPc)):
        return
    L += ["## Sensitivity (A2.5)", "", "✓: false share of the top 1% below liana6's in the same run; ✗: not below.", ""]
    for name, res in (("Calibration", sc), ("Permutations, calibration (receptors, rep 0)", sPc)):
        if res:
            L += [f"{name}:", ""]
            L += table(["Run", *SHOWN], [[t] + [f"{r['false_share_top'][v]:.3f}" + ("" if v == BASELINE else
                                                 " ✓" if r["false_share_top"][v] < r["false_share_top"][BASELINE] else " ✗")
                                                 for v in SHOWN] for t, r in res["runs"].items()])
            n = len(res["runs"])
            L += [f"{name}: runs below liana6 — " + ", ".join(
                f"{v} {sum(r['false_share_top'][v] < r['false_share_top'][BASELINE] for r in res['runs'].values())} of {n}"
                for v in SHOWN if v != BASELINE) + ".", ""]
    if sp:
        L += ["Power, f = 1.5 (power in the top 1% / AUROC):", ""]
        L += table(["Run", *SHOWN], [[t] + [f"{r[v]['single']['power_top']:.2f} / {r[v]['single']['auroc']:.3f}" for v in SHOWN]
                                    for t, r in sp["runs"].items()])
    for name, res in (("Seeds, real labels (Mammary), against seed 1337", ss), ("5,000 vs 1,000 permutations, real labels", sP)):
        if res:
            L += [f"{name} (Spearman of specificity_rank / top-1% Jaccard):", ""]
            L += table(["Run", *SHOWN], [[t] + [f"{r[v]['spearman']:.4f} / {r[v]['top1_jaccard']:.3f}" if v in r else "not run"
                                                for v in SHOWN] for t, r in res["runs"].items()])
    if ss or sP:
        L += ["\"not run\": A2.5's seed and permutation runs compute no MDT (`b2=False`), so the MDT variants have no "
              "second table to compare (as in AMENDMENT_7).", ""]


def cases(L):
    cs = f8("cases")
    if not cs:
        return
    L += ["## Case studies (A2.9): what swap4_<n> drops among liana6's calls (s ≠ t)", "",
          "Breadth: share of cell types where the gene is detected in ≥ 10% of cells; a row's breadth is the larger of its "
          "ligand's and receptor's. p: one-sided Mann–Whitney, dropped broader than kept. Overlap with swap3_<n>'s drops. "
          "Computed locally from this amendment's real tables.", ""]
    rows = []
    for name, r in cs["runs"].items():
        for v in SWAP4[1:]:
            x = r[v]
            rows.append([SHORT[name], v, x["kept"], x["dropped"],
                         f"{x['breadth_dropped_median']:.3f} / {x['breadth_kept_median']:.3f}" if x["breadth_dropped_median"] is not None else "–",
                         "–" if x["mann_whitney_p"] is None else f"{x['mann_whitney_p']:.2g}",
                         f"{x['dropped_by_both']} / {x['dropped_swap4_only']} / {x['dropped_swap3_only']}"])
    L += table(["Data set", "Variant", "kept", "dropped", "breadth dropped / kept", "p", "dropped by both / swap4 only / swap3 only"], rows)


def tool_level(L):
    tl, tp, tr = f8("tool_labels"), f8("tool_power"), f8("tool_real")
    if tl or tp or tr:
        L += ["## Benchmark 1's tests at tool level, inside LIANA+, all four tools", "",
              "p from each null, BH at 5% over the rows with a p. JLP twice: LIANA+'s p (#/B) and Benchmark 1's "
              "((1 + #)/(B + 1)). NATMI's p are on `spec_weight`.", ""]
    if tl:
        L += ["### Labels shuffled (any call is false): runs with any call / mean number of calls", ""]
        rows = []
        for name in dict.fromkeys(r["dataset"] for r in tl["runs"].values()):
            rs = [r for r in tl["runs"].values() if r["dataset"] == name]
            for tool in TOOLS:
                rows.append([SHORT[name], tool, len(rs)] + [f"{sum(r[tool][n]['any_call'] for r in rs)} / "
                                                            f"{np.mean([r[tool][n]['calls'] for r in rs]):,.1f}" for n in TOOL_NULLS])
        L += table(["Data set", "Tool", "reps", *(NAME[n] for n in TOOL_NULLS)], rows)
    if tp:
        L += ["### Power (F2's design, 30 single pairs, f = 1.5, Mammary): power / FDP", ""]
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
        tests = [(k, t, n) for k, r in tr["runs"].items() for t in TOOLS if "spatial" in r.get(t, {})
                 for n in r[t]["spatial"] if r[t]["spatial"][n].get("odds_ratio") is not None]
        q = dict(zip(tests, bh(np.array([tr["runs"][k][t]["spatial"][n]["p"] for k, t, n in tests])))) if tests else {}
        if tests:
            L += ["### Spatial truth among JLP's calls ((1 + #)/(B + 1)): odds ratio kept vs dropped [95% CI]; q: BH over all "
                  "tests below (each session wrote q over its own data set; recomputed here)", ""]
            rows = []
            for k, r in tr["runs"].items():
                for t in TOOLS:
                    if "spatial" in r.get(t, {}):
                        rows.append([SHORT[k], t] + [(lambda x, key: f"{x['odds_ratio']:.2f} [{x['ci95'][0]:.2f}, {x['ci95'][1]:.2f}] "
                                                      f"q {q[key]:.2f}" if x.get("odds_ratio") is not None else "–")(
                            r[t]["spatial"].get(n, {}), (k, t, n)) for n in TOOL_NULLS[2:]])
            L += table(["Data set", "Tool", *(NAME[n] for n in TOOL_NULLS[2:])], rows)
            sig = [f"{SHORT[k]} {t} {NAME[n]}" for (k, t, n), x in q.items() if x < 0.05]
            L += [f"Tests with q < 0.05: {', '.join(sig) if sig else 'none'}.", ""]


def checks(L):
    L += ["## Checks", "",
          "liana6, swap_<n> and swap3_<n> against AMENDMENT_6 (same seeds); CellPhoneDB's, CellChat's and Connectome's "
          "numbers against AMENDMENT_7; NATMI's LOO-t against Connectome's (AMENDMENT_8 expected them equal; see below).", ""]
    fs = lambda r, v: r["false_share_top"].get(v)
    pw = lambda r, v: (r[v]["single"]["power_top"], r[v]["single"]["auroc"]) if v in r else None
    st = lambda r, v: (r[v]["spearman"], r[v]["top1_jaccard"]) if v in r else None
    rows = []
    for n, g, old in (("calibration", fs, OUT6), ("confirm_calibration", fs, OUT6), ("mouse_calibration", fs, OUT6),
                      ("power", pw, OUT6), ("hard_power", pw, OUT6), ("confirm_power", pw, OUT6),
                      ("sensitivity_calibration", fs, OUT7), ("sensitivity_power", pw, OUT7), ("sensitivity_seed", st, OUT7),
                      ("sensitivity_permutations", st, OUT7), ("sensitivity_permutations_calibration", fs, OUT7)):
        rows.append([f"{n} (vs {old.name})", same_runs(f8(n), J(old / f"{n}.json"), g, EARLIER)])
    real, r6 = f8("real"), J(OUT6 / "real.json")
    if real and r6:
        both = [t for t in real["runs"] if t in r6["runs"]]
        n = sum(all(real["runs"][t]["calls"][v] == r6["runs"][t]["calls"][v] for v in EARLIER if v in r6["runs"][t]["calls"]) for t in both)
        rows.append(["real calls (vs amendment6)", f"{n} of {len(both)} data sets identical"])
    for part, keys in (("donor_null", lambda r, t: {k: v for k, v in r[t].items() if k != "rows"}),
                       ("tool_labels", lambda r, t: r[t]), ("tool_power", lambda r, t: r[t]),
                       ("tool_real", lambda r, t: {k: v for k, v in r[t].items() if k != "spatial"})):
        new, old = f8(part), J(OUT7 / f"{part}.json")
        if new and old:
            both = [t for t in new["runs"] if t in old["runs"]]
            n = sum(all(keys(new["runs"][t], tool) == keys(old["runs"][t], tool) for tool in TOOLS[:3]) for t in both)
            rows.append([f"{part}: CellPhoneDB, CellChat, Connectome (vs amendment7)", f"{n} of {len(both)} runs identical"])
    dn, tr = f8("donor_null"), f8("tool_real")
    if dn:
        gap = max(abs(r["NATMI"][f"t2_le_{a}"] - r["Connectome"][f"t2_le_{a}"]) for r in dn["runs"].values() for a in (0.01, 0.05))
        n = sum(all(r["NATMI"][f"t2_le_{a}"] == r["Connectome"][f"t2_le_{a}"] for a in (0.01, 0.05)) for r in dn["runs"].values())
        rows.append(["donor_null: NATMI's LOO-t vs Connectome's", f"identical in {n} of {len(dn['runs'])} runs; share of rows with p ≤ "
                     f"0.01 or 0.05 differs by at most {100 * gap:.4f} percentage points"])
    if tr:
        d = [abs(r["NATMI"]["calls_t2"] - r["Connectome"]["calls_t2"]) for r in tr["runs"].values()]
        rows.append(["tool_real: NATMI's LOO-t calls vs Connectome's", f"identical in {sum(x == 0 for x in d)} of {len(d)} data sets; "
                     f"at most {max(d):,} calls apart"])
    L += table(["Check", "Result"], rows)
    L += ["AMENDMENT_8 stated that NATMI's LOO-t equals Connectome's, since within an interaction the sums over types are "
          "constant. That holds for single-gene interactions only (there the two p-values agree to rounding: LIANA+'s "
          "scores are float32). For a complex, LIANA+ uses the subunit that is lowest in each (s, t), so the sums, and "
          "with them NATMI's LOO-t, change across the interaction's rows. On the 600-cell toy data set the two p-values "
          "differ by at most 6e-8 on single-gene rows and by up to 0.22 on 32 of 80 complex rows.", ""]


def write():
    L = ["# AMENDMENT_8 results: NATMI's specificity replaced too, through every test", "",
         "Plan: [`AMENDMENT_8.md`](../../../liana_final/AMENDMENT_8.md). `swap_<n>`: CellPhoneDB's and CellChat's p replaced "
         "by null n's; `swap3_<n>`: Connectome's too; `swap4_<n>`: NATMI's too. Names as in the paper: JLP = Null A, "
         "BPT = Null B, MDT = B2, LOO-t = T2. Intervals: 95% t-intervals over replicates.", "",
         "Where it ran: 31 of the 42 sessions in `kaggle/kernels.py` on Kaggle; the other 11 (`a8-sensitivity-*`, "
         "`a8-donor-mammary`, `a8-donor-tnbc`, `a8-tool-real-*`) on a lab server, from the same commit (be09388) and "
         "lock file, one CPU each; every run records its machine. Biology and case studies ran locally from the real-data "
         "tables. The tool-level spatial q-values are recomputed here over every data set (each session wrote q over its own).", ""]
    consensus(L)
    donor(L)
    sensitivity(L)
    cases(L)
    tool_level(L)
    checks(L)
    OUT8.mkdir(parents=True, exist_ok=True)
    (OUT8 / "README.md").write_text("\n".join(L) + "\n")
    print("wrote", OUT8 / "README.md")


if __name__ == "__main__":
    write()
