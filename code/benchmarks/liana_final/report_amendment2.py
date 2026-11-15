"""results/liana_final/amendment2/README.md from AMENDMENT_2's JSON files (parts not run yet are skipped)."""
from __future__ import annotations

import json

import numpy as np

from benchmarks.liana_nulls.stages import mean_ci

from .amendment2 import FILTERS, OUT2
from .pipeline import BASELINE, OUT

SHOW = ["liana6", "swap_null_b", "extra_null_b", "swap_b2", "liana6_prop20", "liana6_prop30", "liana6_de",
        "liana6@mean", "swap_null_b@mean"]
LABEL = {"liana5": "LIANA+ as shipped", "liana6": "**liana6 (baseline)**", "swap_null_b": "**swap: Null B**",
         "extra_null_b": "**7th/8th: Null B**", "swap_b2": "swap: B2", "extra_b2": "7th/8th: B2", "extra_null_a": "7th/8th: Null A",
         "liana6_prop20": "filter: detected ≥ 20%", "liana6_prop30": "filter: detected ≥ 30%", "liana6_de": "filter: both sides DE",
         "liana6@mean": "liana6, mean aggregation", "swap_null_b@mean": "swap: Null B, mean aggregation",
         "extra_null_b@mean": "7th/8th: Null B, mean aggregation", "swap_null_b_donor": "swap: Null B, donor permutations"}
DS = {"tnbc_wu2021": "TNBC", "allen_brain_atlas": "Allen brain", "tabula_sapiens_mammary": "TS Mammary",
      "tabula_sapiens_liver": "TS Liver", "tabula_sapiens_lung": "TS Lung", "moffitt2018_hypothalamus": "Hypothalamus"}


def J(name):
    f = OUT2 / f"{name}.json"
    return json.loads(f.read_text()) if f.exists() else None


def fm(m, d=3):
    v, (lo, hi) = m
    return f"{v:.{d}f}" if lo == hi else f"{v:.{d}f} [{lo:.{d}f}, {hi:.{d}f}]"


def cal_summary(runs: dict) -> dict:
    """{(dataset, setting, scenario): {variant: mean_ci, "all": mean_ci}}"""
    groups = {}
    for r in runs.values():
        groups.setdefault((r["dataset"], r.get("setting", "default"), r["scenario"]), []).append(r)
    return {k: {"all": mean_ci([r["false_share_all"] for r in rs]), "n": len(rs)}
            | {v: mean_ci([r["false_share_top"][v] for r in rs]) for v in rs[0]["false_share_top"]}
            for k, rs in groups.items()}


def cal_table(S: dict, variants=SHOW) -> list:
    keys = sorted(S)
    L = ["| | " + " | ".join(f"{DS[d]}{'' if s == 'default' else ' ' + s}, {sc}s shuffled (n={S[(d, s, sc)]['n']})"
                             for d, s, sc in keys) + " |", "|---|" + "---|" * len(keys),
         "| all rows | " + " | ".join(f"{S[k]['all'][0]:.3f}" for k in keys) + " |"]
    for v in variants:
        if all(v in S[k] for k in keys):
            L.append(f"| {LABEL.get(v, v)} | " + " | ".join(fm(S[k][v]) for k in keys) + " |")
    return L


def power_summary(runs: dict) -> dict:
    """{(dataset, setting, design, f): {variant: {class: {power_top, auroc}}, precision_at_k}}"""
    groups = {}
    for r in runs.values():
        groups.setdefault((r["dataset"], r.get("setting", "default"), r["design"], r["effect"]), []).append(r)
    out = {}
    for k, rs in groups.items():
        vs = [v for v in rs[0] if isinstance(rs[0][v], dict) and "precision_at_k" in rs[0][v]]
        classes = sorted(rs[0]["planted_rows_present"])
        out[k] = {"n": len(rs), "classes": classes,
                  "present": {c: float(np.mean([r["planted_rows_present"][c] for r in rs])) for c in classes}}
        for v in vs:
            out[k][v] = {"precision_at_k": mean_ci([r[v]["precision_at_k"] for r in rs])}
            for c in classes:
                out[k][v][c] = {m: mean_ci([r[v][c][m] for r in rs]) for m in ("power_top", "auroc")}
            if v != BASELINE:
                out[k][v]["gap_single"] = mean_ci([r[v]["single"]["power_top"] - r[BASELINE]["single"]["power_top"] for r in rs])
    return out


def power_table(P: dict, keys: list, cls: str, variants=SHOW) -> list:
    L = ["| | " + " | ".join(f"{DS[d]}{'' if s == 'default' else ' ' + s}, f = {f} (n={P[(d, s, de, f)]['n']})"
                             for d, s, de, f in keys) + " |", "|---|" + "---|" * len(keys)]
    for v in variants:
        if all(v in P[k] for k in keys):
            L.append(f"| {LABEL.get(v, v)} | " + " | ".join(
                f"{P[k][v][cls]['power_top'][0]:.2f} / {P[k][v][cls]['auroc'][0]:.3f} / {P[k][v]['precision_at_k'][0]:.2f}"
                for k in keys) + " |")
    return L


def confirm_verdict(S: dict, P: dict) -> dict:
    """AMENDMENT_1's criteria 1 and 2 on Lung."""
    out = {}
    for v in ("swap_null_b", "extra_null_b"):
        cal = [S[k][v][0] < S[k][BASELINE][0] for k in S if k[0] == "tabula_sapiens_lung" and k[1] == "default"]
        pk = [k for k in P if k[0] == "tabula_sapiens_lung" and k[3] == 2.0]
        c1 = all(cal) if len(cal) == 2 else None
        c2 = None
        if pk:
            m, b = P[pk[0]][v]["single"], P[pk[0]][BASELINE]["single"]
            c2 = m["power_top"][0] >= 0.8 and m["auroc"][0] >= b["auroc"][0]
        out[v] = dict(calibration=c1, power=c2, verdict=None if None in (c1, c2) else
                      ("confirmed" if c1 and c2 else "not confirmed (" + ", ".join(n for n, c in (("calibration", c1), ("power", c2)) if not c) + ")"))
    return out


def write():
    L = ["# AMENDMENT_2: robustness, held-out confirmation and reproducibility", "",
         "Pre-registered in `benchmarks/liana_final/AMENDMENT_2.md` (committed before its runs). Run: "
         "`python -m benchmarks.liana_final.amendment2 --part all`. Measures as in AMENDMENT_1: false share in the top 1% "
         "(lower is better); power in the top 1% / AUROC planted vs one-sided / precision at k. Means over replicates "
         "[95% t-interval].", ""]
    yn = lambda x: "pending" if x is None else ("yes" if x else "no")
    # A2.1
    cc, cp = J("confirm_calibration"), J("confirm_power")
    if cc or cp:
        L += ["## A2.1 Held-out confirmation: Tabula Sapiens Lung", ""]
        S = cal_summary(cc["runs"]) if cc else {}
        P = power_summary(cp["runs"]) if cp else {}
        if S:
            L += ["Calibration:", ""] + cal_table(S) + [""]
        if P:
            L += ["Power (single pairs): power in the top 1% / AUROC / precision at k:", ""] + power_table(P, sorted(P), "single") + [""]
        V = confirm_verdict(S, P)
        L += ["Frozen rule (AMENDMENT_1, criteria 1 and 2):", "", "| Variant | Calibration | Power | Verdict |", "|---|---|---|---|"]
        L += [f"| {LABEL[v]} | {yn(d['calibration'])} | {yn(d['power'])} | {d['verdict'] or 'pending'} |" for v, d in V.items()] + [""]
    # A2.2
    hp = J("hard_power")
    if hp:
        P = power_summary(hp["runs"])
        keys = sorted(P)
        L += ["## A2.2 Harder power (TS Mammary background)", "",
              f"Planted rows present: " + "; ".join(f"f = {k[3]}: " + ", ".join(f"{c} {P[k]['present'][c]:.1f}" for c in P[k]["classes"])
                                                    for k in keys) + ".", ""]
        for cls, n in (("single", 30), ("complex", 10), ("unequal", 10)):
            L += [f"**{cls}** ({n} pairs a replicate): power in the top 1% / AUROC planted vs one-sided / precision at k (all true rows):", ""]
            L += power_table(P, keys, cls) + [""]
        L += ["Power gap to liana6 (single pairs, power in the top 1%):", "",
              "| | " + " | ".join(f"f = {k[3]}" for k in keys) + " |", "|---|" + "---|" * len(keys)]
        for v in ("swap_null_b", "extra_null_b", "liana6_de"):
            L.append(f"| {LABEL[v]} | " + " | ".join(fm(P[k][v]["gap_single"], 2) for k in keys) + " |")
        L.append("")
    # A2.3
    dn, dr = J("donor_null"), J("donor_real")
    if dn:
        L += ["## A2.3 Donor-aware permutations", "", "(a) Labels shuffled within donor: share of rows (s ≠ t, score > 0) "
              "with p ≤ 0.05 (nominal 0.05), mean over replicates.", "",
              "| Dataset (donors) | Tool | Null A, global | Null A, within donor | Null B, global | Null B, within donor |",
              "|---|---|---|---|---|---|"]
        runs = dn["runs"].values()
        for d in dict.fromkeys(r["dataset"] for r in runs):
            for tool in ("CellPhoneDB", "CellChat"):
                g = {s: [r for r in runs if r["dataset"] == d and r["scheme"] == s] for s in ("global", "donor")}
                if not g["global"] or not g["donor"]:
                    continue
                c = lambda s, n: fm(mean_ci([r[tool][f"{n}_le_0.05"] for r in g[s]]))
                L.append(f"| {DS[d]} ({g['global'][0]['n_donors']}) | {tool} | {c('global', 'null_a')} | {c('donor', 'null_a')} | "
                         f"{c('global', 'null_b')} | {c('donor', 'null_b')} |")
        L.append("")
    if dr:
        L += ["(b) Real labels: swap_null_b with within-donor against global permutations.", "",
              "| Dataset | Top-1% Jaccard | Calls, global | Calls, within donor | Call Jaccard |", "|---|---|---|---|---|"]
        for d, r in dr["runs"].items():
            x = r["swap_null_b"]
            L.append(f"| {DS[d]} | {x['top1_jaccard']:.3f} | {x['calls_global']:,} | {x['calls_donor']:,} | {x['call_jaccard']:.3f} |")
        L.append("")
    # A2.4 head to head
    allcal = {k: v for n in ("confirm_calibration", "recalibration", "mouse_calibration", "sensitivity_calibration")
              if (j := J(n)) for k, v in cal_summary(j["runs"]).items()}
    allpow = {k: v for n in ("confirm_power", "hard_power", "sensitivity_power") if (j := J(n)) for k, v in power_summary(j["runs"]).items()}
    if allcal or allpow:
        L += ["## A2.4 Null B against simple filters", "",
              "Settings (calibration scenarios and power settings of this amendment) in which the Null B variant beats the filter: "
              "lower false share in the top 1% / higher AUROC (single pairs).", "",
              "| | " + " | ".join(LABEL[f] for f in FILTERS) + " |", "|---|" + "---|" * len(FILTERS)]
        for v in ("swap_null_b", "extra_null_b"):
            cells = []
            for f in FILTERS:
                c = [allcal[k][v][0] < allcal[k][f][0] for k in allcal if f in allcal[k] and v in allcal[k]]
                p = [allpow[k][v]["single"]["auroc"][0] > allpow[k][f]["single"]["auroc"][0] for k in allpow if f in allpow[k]]
                cells.append(f"calibration {sum(c)}/{len(c)}; AUROC {sum(p)}/{len(p)}")
            L.append(f"| {LABEL[v]} | " + " | ".join(cells) + " |")
        L.append("")
    # A2.5
    sc, sp_, ss, sperm, spc = (J(n) for n in ("sensitivity_calibration", "sensitivity_power", "sensitivity_seed",
                                              "sensitivity_permutations", "sensitivity_permutations_calibration"))
    if any((sc, sp_, ss, sperm, spc)):
        L += ["## A2.5 Sensitivity (TS Mammary)", ""]
        if ss:
            L += ["Seed (real labels, against seed 1337): Spearman of specificity_rank / Jaccard of the top 1%.", "",
                  "| Seed | liana6 | swap: Null B | 7th/8th: Null B |", "|---|---|---|---|"]
            L += [f"| {t[4:]} | " + " | ".join(f"{r[v]['spearman']:.4f} / {r[v]['top1_jaccard']:.3f}" for v in ("liana6", "swap_null_b", "extra_null_b"))
                  + " |" for t, r in ss["runs"].items()] + [""]
        if sperm or spc:
            L += ["Permutations, 1,000 vs 5,000 (the hypothalamus, and Mammary):", ""]
            if sperm and "real" in sperm["runs"]:
                r = sperm["runs"]["real"]
                L += ["| Real labels | liana6 | swap: Null B | 7th/8th: Null B |", "|---|---|---|---|",
                      "| Spearman / top-1% Jaccard | " + " | ".join(f"{r[v]['spearman']:.4f} / {r[v]['top1_jaccard']:.3f}"
                                                         for v in ("liana6", "swap_null_b", "extra_null_b")) + " |", ""]
            if spc:
                L += cal_table(cal_summary(spc["runs"]), ["liana6", "swap_null_b", "extra_null_b"]) + [""]
        if sc:
            L += ["Calibration under each setting (rep 0):", ""] + cal_table(cal_summary(sc["runs"])) + [""]
        if sp_:
            P = power_summary(sp_["runs"])
            L += ["Power under each setting (f = 1.5, single pairs):", ""] + power_table(P, sorted(P), "single") + [""]
    # A2.6
    rc, old = J("recalibration"), (json.loads((OUT / "calibration.json").read_text()) if (OUT / "calibration.json").exists() else None)
    if rc:
        L += ["## A2.6 The p-values themselves", "", "Calibration rerun with the same seeds (Mammary, Liver):", ""]
        L += cal_table(cal_summary(rc["runs"])) + [""]
        if old:
            same, n = 0, 0
            for t, r in rc["runs"].items():
                d, scen, rep = t.split("|")
                o = old["results"].get(d, {}).get(scen, {}).get("runs", [])
                if int(rep) < len(o):
                    n += 1
                    same += all(abs(o[int(rep)][v]["false_share_top"] - r["false_share_top"][v]) < 1e-12 for v in o[int(rep)] if v in r["false_share_top"])
            L += [f"Reproducibility: {same} of {n} reruns give the false shares of `calibration.json` exactly.", ""]
        L += ["Share of false rows (shuffled side) with p ≤ α, pooled over runs; BH at 5% within each tool: share of discoveries that "
              "are false rows (beside the share of false rows overall).", "",
              "| Tool | Null | p ≤ 0.01 | p ≤ 0.05 | p ≤ 0.1 | BH discoveries | false share of discoveries | false share of all rows |",
              "|---|---|---|---|---|---|---|---|"]
        runs = [r for r in rc["runs"].values() if "p_calibration" in r]
        for tool in ("CellPhoneDB", "CellChat"):
            for n in ("null_a", "null_b", "b2"):
                at = {a: np.mean([r["p_calibration"][tool][n]["at_most"][a] for r in runs]) for a in ("0.01", "0.05", "0.1")}
                disc = np.mean([r["p_calibration"][tool][n]["bh_discoveries"] for r in runs])
                fs = [r["p_calibration"][tool][n]["bh_false_share"] for r in runs if r["p_calibration"][tool][n]["bh_false_share"] is not None]
                fa = np.mean([r["p_calibration"][tool]["false_share_all"] for r in runs])
                nm = {"null_a": "Null A", "null_b": "Null B", "b2": "B2"}[n]
                L.append(f"| {tool} | {nm} | {at['0.01']:.4f} | {at['0.05']:.4f} | {at['0.1']:.4f} | {disc:,.0f} | "
                         f"{np.mean(fs) if fs else float('nan'):.3f} | {fa:.3f} |")
        L.append("")
    # A2.7
    mc = J("mouse_calibration")
    if mc:
        S = cal_summary(mc["runs"])
        L += ["## A2.7 A second species: Allen brain (mouse)", ""] + cal_table(S) + [""]
        L += ["AMENDMENT_1's criterion 1 here: " + ", ".join(
            f"{LABEL[v]} {'yes' if all(S[k][v][0] < S[k][BASELINE][0] for k in S) else 'no'}" for v in ("swap_null_b", "extra_null_b")), ""]
    # A2.8
    bio = J("biology")
    if bio and "tests" in bio:
        T = bio["tests"]
        L += ["## A2.8 Biology: every test, multiplicity and detectable effects", "",
              "Odds ratio among the baseline's calls, kept vs dropped [95% cluster-bootstrap interval]; two-sided bootstrap p; BH q within "
              "the variant and across all; detectable: the odds ratios each test could detect with 80% power (α = 0.05).", ""]
        for v in ("swap_null_b", "extra_null_b", "swap_null_b_donor", "swap_b2", "extra_b2", "extra_null_a", "liana5"):
            rows = [(k.split(" | ")[1], r) for k, r in T.items() if k.split(" | ")[0] == v]
            if not rows:
                continue
            L += [f"**{LABEL[v]}**", "", "| Test | OR [95% CI] | p | q (variant) | q (all) | detectable |", "|---|---|---|---|---|---|"]
            for t, r in rows:
                if r["odds_ratio"] is None:
                    L.append(f"| {t} | – (keeps or drops none) | | | | |")
                else:
                    lo, hi = r["ci95"]
                    L.append(f"| {t} | {r['odds_ratio']:.2f} [{lo:.2f}, {hi:.2f}] | {r['p']:.4f} | {r.get('q_variant', float('nan')):.3f} | "
                             f"{r.get('q_all', float('nan')):.3f} | ≤ {r['detectable'][0]:.2f} or ≥ {r['detectable'][1]:.2f} |")
            L.append("")
    # A2.9
    cs = J("cases")
    if cs:
        L += ["## A2.9 What Null B drops, and why", "",
              "Among liana6's calls (s ≠ t): rows swap_null_b drops vs keeps. Limiting side: the side with the larger Null B p. Breadth: "
              "share of cell types in which the limiting side's gene is detected in ≥ 10% of cells. Mann–Whitney: dropped broader than kept.", "",
              "| Dataset | Tool | Dropped | Kept | Limiting side of dropped (ligand / receptor) | Breadth, dropped vs kept (median) | p |",
              "|---|---|---|---|---|---|---|"]
        for d, r in cs["runs"].items():
            for tool in ("CellPhoneDB", "CellChat"):
                x = r[tool]
                L.append(f"| {DS[d]} | {tool} | {x['dropped']:,} | {x['kept']:,} | {x['limiting_side_dropped']['ligand']:,} / "
                         f"{x['limiting_side_dropped']['receptor']:,} | {x['breadth_dropped_median']:.2f} vs {x['breadth_kept_median']:.2f} | "
                         f"{x['mann_whitney_p']:.2g} |")
        L += ["", "Top dropped and promoted rows: `cases_<dataset>_<tool>_dropped.csv`, `cases_<dataset>_promoted.csv`.", ""]
    # A2.10
    tm = J("timing")
    if tm:
        L += ["## A2.10 Run time", "", "| Dataset | Cells | li.mt.rank_aggregate (s) | This pipeline (s) | Ratio |", "|---|---|---|---|---|"]
        L += [f"| {DS[d]} | {r['n_cells']:,} | {r['liana_rank_aggregate']:,.0f} | {r['this_pipeline']:,.0f} | {r['ratio']:.2f} |"
              for d, r in tm["runs"].items()] + ["", "This pipeline: six methods (LIANA+'s five and CellChat), Null B, and B2 with its decoy genes.", ""]
    def machines(parts, keep=lambda tag: True):
        return sorted({r["machine"] for p in parts if J(p) for t, r in J(p)["runs"].items() if "machine" in r and keep(t)})
    lung = machines(("confirm_calibration", "confirm_power", "donor_null"), lambda t: t.startswith("tabula_sapiens_lung"))
    perm = machines(("sensitivity_permutations_calibration",), lambda t: "tabula_sapiens_mammary@n_perms=5000" in t)
    L += ["## Deviations from AMENDMENT_2", "",
          "- **Machines.** AMENDMENT_2 was written for one 8 GB laptop. Every result here was produced on Kaggle (4 cores, "
          "32 GB) or, where noted, a lab server, from one commit in the same locked environment "
          "(`requirements.txt`), with the same code, seeds and settings. Each run records its machine; "
          "`kaggle/README.md` lists the sessions and `RERUNS.md` compares them with earlier laptop runs."]
    if lung:
        L += ["- **Lung (A2.1, A2.3): all cells, not a subsample.** A full-Lung run was killed by the operating system for "
              "memory on the laptop, and so was the pre-registered fallback of 30,000 cells: memory is set by the permutation "
              "tensor, 2 × 1,000 × rows, and rows grow with the square of the number of cell types (34 in Lung, about 370,000 "
              "rows), not with the number of cells. Every Lung run used all 65,847 cells, on " + ", ".join(lung) + "."]
    L += ["- **Permutations (A2.5).** 5,000 permutations on Tabula Sapiens Mammary need a 2 × 5,000 × 86,000 permutation "
          "tensor (about 7 GB) per tool, beyond the laptop. "
          + ("The Mammary comparison (real labels, receptor-shuffled calibration rep 0) ran on a lab server with more "
             "memory than a Kaggle session (" + ", ".join(perm) + "); the hypothalamus on Kaggle." if perm else
             "The comparison was run on the hypothalamus (10,197 rows), real labels and the receptor-shuffled calibration (rep 0).")]
    L += [""]
    OUT2.mkdir(parents=True, exist_ok=True)
    (OUT2 / "README.md").write_text("\n".join(L) + "\n")
    print("wrote", OUT2 / "README.md")
