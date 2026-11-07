"""results/liana_nulls/README.md from the JSON files, with the plan's decision rule applied as written."""
from __future__ import annotations

import json

from .common import OUT, amendment_hash, plan_hash

DATASETS = ["tnbc_wu2021", "allen_brain_atlas", "tabula_sapiens_mammary", "tabula_sapiens_liver", "moffitt2018_hypothalamus"]
TOOLS = ["cellchat", "cellphonedb"]


def _load(name):
    f = OUT / f"{name}.json"
    return json.loads(f.read_text()) if f.exists() else None


def _ci(x, d=3):
    return f"[{x[0]:.{d}f}, {x[1]:.{d}f}]"


def fmt(v):
    return "–" if v is None else (f"{v:,}" if isinstance(v, int) else str(v))


def fidelity_and_columns(lines):
    lines += ["## T0 fidelity and the three columns", "",
              "Column 1 is the tool itself (CellChat 2.1.2 in R; CellPhoneDB 5.0.1). Column 2 is our Python Null A on "
              "the tool's own permutations. Column 3 is Null B. Calls: the tool's own rule for columns 1 and 2, "
              "BH q ≤ 0.05 over the same rows for Null A and Null B.", "",
              "| Dataset | Tool | Rows where p (col. 1) = p (col. 2) | Calls, tool's rule: col. 1 / col. 2 | Null A (BH) | "
              "Null B (BH) | Null B outside Null A | Null A only: ligand / receptor fails | LIANA+ built-in vs tool rule (Jaccard) |",
              "|---|---|---|---|---|---|---|---|---|"]
    for ds in DATASETS:
        c = _load(f"core_{ds}")
        if not c:
            continue
        for t in TOOLS:
            s = c.get(t)
            if not s:
                continue
            fd, k = s.get("fidelity") or {}, s["calls"]
            same = f"{fd.get('pval_identical', 0):,} / {fd.get('n_rows', 0):,}" if fd else "–"
            native = fd.get("calls_r", fd.get("calls_native")) if fd else None
            port = fd.get("calls_python") if fd else None
            lo = s["null_a_only_by_limiting_partner"]
            lb = s.get("liana_builtin", {})
            lines.append(f"| {ds}{' (Amendment 1)' if ds == 'moffitt2018_hypothalamus' else ''} | {t} | {same} | {fmt(native)} / {fmt(port)} | {k['2_our_null_a_bh']:,} | {k['3_null_b_bh']:,} | "
                         f"{s['null_b_calls_outside_null_a']} | {lo.get('ligand', 0):,} / {lo.get('receptor', 0):,} | "
                         f"{lb.get('jaccard_with_native_rule', float('nan')):.3f} |")
    lines.append("")


def calibration(lines):
    c = _load("calibration")
    if not c:
        return {}
    lines += ["## T1 calibration (Tabula Sapiens, known-null data, α = 0.05)", "",
              "| Tool | Tissue | Receptor genes shuffled: FDR Null A / Null B | Ligand genes shuffled: FDR Null A / Null B | "
              "Shuffled labels: P(any call) A / B |", "|---|---|---|---|---|"]
    ok = True
    for t, by in c["results"].items():
        for tissue, s in by.items():
            r, l, lab = s["receptor"], s["ligand"], s["labels"]
            lines.append(f"| {t} | {tissue} | {r['null_a']['fdr']:.3f} / **{r['null_b']['fdr']:.3f}** {_ci(r['null_b']['ci95'])} | "
                         f"{l['null_a']['fdr']:.3f} / **{l['null_b']['fdr']:.3f}** {_ci(l['null_b']['ci95'])} | "
                         f"{lab['null_a']['p_any_call']:.2f} / {lab['null_b']['p_any_call']:.2f} |")
            ok &= r["null_b"]["ci95"][1] <= 0.05 and l["null_b"]["ci95"][1] <= 0.05
    lines.append("")
    return {"t1_null_b_upper_ci_le_0.05": ok}


def power(lines):
    p = _load("power")
    if not p:
        return {}
    lines += ["## T2 power on planted edges (Tabula Sapiens Mammary, labels shuffled, 30 planted pairs)", "",
              "| Tool | Fold change | Tool's rule: power / FDR | Null A: power / FDR | Null B: power / FDR | "
              "One-sided false calls A / B | Null B calls (BH) |", "|---|---|---|---|---|---|---|"]
    ok, vacuous = True, []
    for t, by in p["summary"].items():
        for f, s in by.items():
            n, a, b = s["native_rule"], s["null_a"], s["null_b"]
            lines.append(f"| {t} | {f} | {n['power'][0]:.2f} / {n['fdr'][0]:.3f} | {a['power'][0]:.2f} / {a['fdr'][0]:.3f} "
                         f"{_ci(a['fdr'][1])} | {b['power'][0]:.2f} / {b['fdr'][0]:.3f} {_ci(b['fdr'][1])} | "
                         f"{a['one_sided_false_calls']:.1f} / {b['one_sided_false_calls']:.1f} | {b['calls']:.1f} |")
            ok &= b["fdr"][1][1] <= 0.05
            if b["calls"] == 0:
                vacuous.append(f"{t} f={f}")
    lines += ["", f"Replicates per cell: {p['design']['reps']}; B = {p['design']['B']}.", ""]
    if vacuous:
        lines += [f"**Vacuous (Amendment 2):** Null B (BH) made no calls at {', '.join(vacuous)}; its FDR there is 0 by "
                  "construction (with B = 1,000 the p-value floor is 1/1001, and BH needs about m/50 rows at the floor).", ""]
    if any("null_b_tool_rule" in s for by in p["summary"].values() for s in by.values()):
        lines += ["### T2 under each tool's own unadjusted rule (Amendment 2)", "",
                  "| Tool | Fold change | Null A, tool's rule: power / FDR [95% CI] | Null B, tool's rule: power / FDR [95% CI] | "
                  "One-sided false calls A / B |", "|---|---|---|---|---|"]
        for t, by in p["summary"].items():
            for f, s in by.items():
                a, b = s["null_a_tool_rule"], s["null_b_tool_rule"]
                lines.append(f"| {t} | {f} | {a['power'][0]:.2f} / {a['fdr'][0]:.3f} {_ci(a['fdr'][1])} | "
                             f"{b['power'][0]:.2f} / {b['fdr'][0]:.3f} {_ci(b['fdr'][1])} | "
                             f"{a['one_sided_false_calls']:.1f} / {b['one_sided_false_calls']:.1f} |")
        lines.append("")
    return {"t2_null_b_upper_ci_le_0.05": ok, "t2_vacuous_at": vacuous}


def biology(lines):
    b = _load("biology")
    if not b:
        return {}
    lines += ["## T3 biology", "",
              "Primary: among Null A's calls (s ≠ t), the odds that a call Null B keeps sits on a colocalised (s, t) "
              "pair, against a call it drops; 95% CI from 2,000 bootstrap resamples of (s, t) pairs.", "",
              "| Dataset | Tool | Kept: colocalised / calls | Dropped: colocalised / calls | Odds ratio [95% CI] |",
              "|---|---|---|---|---|"]
    above, below, worse = False, False, []
    for ds, by in b.items():
        for t in TOOLS:
            if t not in by:
                continue
            k = by[t]["primary_kept_vs_dropped"]
            lines.append(f"| {ds} | {t} | {k['kept_colocalised']:,} / {k['kept_calls']:,} ({k['frac_colocalised_kept']:.1%}) | "
                         f"{k['dropped_colocalised']:,} / {k['dropped_calls']:,} ({k['frac_colocalised_dropped']:.1%}) | "
                         f"{k['odds_ratio']:.2f} {_ci(k['ci95'], 2)} |")
            above |= k["ci95"][0] > 1
            below |= k["ci95"][1] < 1
    lines += ["", "Secondary (AUPRC; Null B − other, paired bootstrap 95% CI):", "",
              "| Dataset | Tool | Truth | Null A | Null B | Tool's rule | Null B − Null A |", "|---|---|---|---|---|---|---|"]
    for ds, by in b.items():
        for t in TOOLS:
            if t not in by:
                continue
            for metric in ("pair_level_auprc", "openproblems_source_target", "cytokine_ligand_target"):
                if metric not in by[t]:
                    continue
                m = by[t][metric]["methods"]
                d = m["null_a"]
                lines.append(f"| {ds} | {t} | {metric} | {m['null_a']['auprc']:.3f} | {m['null_b']['auprc']:.3f} | "
                             f"{m['native_rule']['auprc']:.3f} | {d['auprc_diff_ref_minus_this']:+.3f} {_ci(d['auprc_diff_ci95'])} |")
                if d["auprc_diff_ci95"][1] < -0.02:
                    worse.append(f"{ds}/{t}/{metric}")
    lines.append("")
    return {"t3_primary_ci_above_1_somewhere": above, "t3_primary_ci_below_1_nowhere": not below,
            "t3_no_secondary_significantly_worse": not worse, "worse": worse}


def _or(d):
    return "–" if not d or d.get("odds_ratio") is None else f"{d['odds_ratio']:.2f} {_ci(d['ci95'], 2)}"


def amendment(lines, t3_above, t3_below):
    """Amendment 1: A1 (hypothalamus spatial truth) and A2 (downstream truths), and the extended verdict."""
    a = _load("amendment_1")
    if not a:
        return {}
    lines += ["## Amendment 1 (added after the plan; `AMENDMENT_1.md`, SHA-256 "
              f"`{amendment_hash()[:16]}…`)", "",
              f"Amendment 2 (`AMENDMENT_2.md`, SHA-256 `{amendment_hash(2)[:16]}…`) adds the T2 tool-rule table above. "
              "Not part of the pre-registered decision above. The amendment was written before any T3 result of this "
              "pipeline; it was not changed after a smoke test of A2 on the Allen brain.", ""]
    a1 = a.get("A1_hypothalamus_spatial", {})
    above, below = t3_above, t3_below
    if a1:
        st = a1.get("spatial_truth", {})
        lines += [f"### A1. Hypothalamus (Moffitt 2018): MERFISH neighbourhood truth, {st.get('colocalised', 0)} of "
                  f"{st.get('pairs', 0)} ordered class pairs colocalised", "",
                  "| Tool | Kept: colocalised / calls | Dropped: colocalised / calls | Odds ratio [95% CI] |", "|---|---|---|---|"]
        for t in TOOLS:
            if t in a1:
                k = a1[t]["primary_kept_vs_dropped"]
                lines.append(f"| {t} | {k['kept_colocalised']:,} / {k['kept_calls']:,} ({k['frac_colocalised_kept']:.1%}) | "
                             f"{k['dropped_colocalised']:,} / {k['dropped_calls']:,} ({k['frac_colocalised_dropped']:.1%}) | "
                             f"{k['odds_ratio']:.2f} {_ci(k['ci95'], 2)} |")
                above |= k["ci95"][0] > 1
                below |= k["ci95"][1] < 1
        lines.append("")
    a2 = a.get("A2_downstream", {})
    d_above, d_below, against = False, False, []
    if a2:
        lines += ["### A2. Downstream response of the receiver (Null A's calls, s ≠ t; responds = z ≥ 1.645)", "",
                  "D1: odds of response, calls Null B keeps vs drops. D2: kept vs dropped because the receptor side fails "
                  "(predicted > 1) or the ligand side fails (no prediction). Bootstrap over (footprint, receiver).", "",
                  "| Dataset | Tool | Truth | Calls in scope | Responding: kept / dropped | D1 | D1 (z ≥ 1) | D2 receptor | D2 ligand |",
                  "|---|---|---|---|---|---|---|---|---|"]
        for ds, by in a2.items():
            for t in TOOLS:
                for res, e in by.get(t, {}).items():
                    d1 = e["D1"]
                    frac = (f"{d1['frac_responding_kept']:.1%} / {d1['frac_responding_dropped']:.1%}"
                            if d1.get("odds_ratio") is not None else "–")
                    lines.append(f"| {ds} | {t} | {res} | {e['n_lines']:,} | {frac} | {_or(d1)} | {_or(e['D1_z1'])} | "
                                 f"{_or(e['D2_receptor'])} | {_or(e['D2_ligand'])} |")
                    if d1.get("ci95"):
                        d_above |= d1["ci95"][0] > 1
                    for k in ("D1", "D2_receptor", "D2_ligand"):     # AMENDMENT_1: "any D1 or D2 CI entirely below 1"
                        if (e[k].get("ci95") or [0, 2])[1] < 1:
                            d_below = True
                            against.append(f"{ds}/{t}/{res}/{k} {_or(e[k])}")
        lines.append("")
    downstream = "against" if d_below else ("supports" if d_above else "neutral")
    lines += ["### Extended verdict (Amendment 1)", "",
              f"- Criterion (2) over all three spatial truths: CI above 1 somewhere = {above}, below 1 somewhere = {below} "
              f"→ **{'met' if above and not below else 'not met'}**.",
              f"- Downstream truths: **{downstream}**" + (f", from: {'; '.join(against)}." if against else "."), ""]
    return {"amendment_spatial_met": above and not below, "amendment_downstream": downstream}


def status(lines):
    """One row per component: pre-registered or an amendment, and whether its result files exist."""
    rows = [("T0 fidelity + columns 1-3", "pre-registered", [f"core_{d}" for d in DATASETS[:4]]),
            ("T1 calibration", "pre-registered", ["calibration"]),
            ("T2 power", "pre-registered", ["power"]),
            ("T2 under the tools' own rule", "Amendment 2", ["power"]),
            ("T3 biology (TNBC, Allen brain)", "pre-registered", ["biology"]),
            ("A1 hypothalamus: core", "Amendment 1", ["core_moffitt2018_hypothalamus"]),
            ("A1 hypothalamus: spatial truth", "Amendment 1", ["amendment_1"]),
            ("A2 downstream truths (CytoSig, PROGENy)", "Amendment 1", ["amendment_1"])]
    from .report3 import STATUS_ROWS
    rows += STATUS_ROWS
    lines += ["## Status", "", "| Component | Plan | Result files | Status |", "|---|---|---|---|"]
    for name, plan, files in rows:
        have = [f for f in files if (OUT / f"{f}.json").exists()]
        if name.startswith("A1 hypothalamus: spatial"):
            ok = bool((_load("amendment_1") or {}).get("A1_hypothalamus_spatial", {}).get("cellchat"))
        else:
            ok = len(have) == len(files)
        lines.append(f"| {name} | {plan} | {len(have)}/{len(files)} | {'done' if ok else 'pending'} |")
    lines.append("")


def write():
    lines = ["# Null A vs Null B in the LIANA+ pipeline: results", "",
             f"Plan: `benchmarks/liana_nulls/PLAN.md` (SHA-256 `{plan_hash()[:16]}…`, fixed before any result). "
             "Run: `python -m benchmarks.liana_nulls`.", ""]
    fidelity_and_columns(lines)
    verdict = {}
    verdict.update(calibration(lines))
    verdict.update(power(lines))
    verdict.update(biology(lines))
    valid = verdict.get("t1_null_b_upper_ci_le_0.05") and verdict.get("t2_null_b_upper_ci_le_0.05")
    if not {"t1_null_b_upper_ci_le_0.05", "t2_null_b_upper_ci_le_0.05", "t3_no_secondary_significantly_worse"} <= set(verdict):
        decision = "incomplete: not every stage has run"
    elif valid and verdict["t3_primary_ci_above_1_somewhere"] and verdict["t3_primary_ci_below_1_nowhere"] \
            and verdict["t3_no_secondary_significantly_worse"]:
        decision = "Null B is the better null for the s → t claim"
    elif valid and verdict["t3_no_secondary_significantly_worse"]:
        decision = "Null B is valid, but not biologically distinguishable from Null A on these truths"
    else:
        decision = "Null B rejected by the plan's rule"
    lines += ["## Decision (the plan's rule, applied as written)", "", f"**{decision}.**", "",
              "```", json.dumps(verdict, indent=1), "```", ""]
    amendment(lines, verdict.get("t3_primary_ci_above_1_somewhere", False), not verdict.get("t3_primary_ci_below_1_nowhere", True))
    if verdict.get("t2_vacuous_at"):
        lines += [f"Note: the T2 part of criterion (1) holds only vacuously at {', '.join(verdict['t2_vacuous_at'])} "
                  "(Amendment 2).", ""]
    from .report3 import amendment3
    amendment3(lines)
    status(lines)
    (OUT / "README.md").write_text("\n".join(lines))
    print("\n".join(lines))
