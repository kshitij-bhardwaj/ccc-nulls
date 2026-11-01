"""The Amendment 3 section of results/liana_nulls/README.md (AMENDMENT_3.md)."""
from __future__ import annotations

from .common import amendment_hash
from .report import DATASETS, TOOLS, _ci, _load

LABEL = {"null_a": "Null A (BH)", "null_b": "**Null B (BH), full**", "null_b|ligand_only": "Null B, ligand side only",
         "null_b|receptor_only": "Null B, receptor side only", "null_b|by": "Null B, BY", "null_b|tool_rule": "Null B, tool's rule",
         "null_a|tool_rule": "Null A, tool's rule", "null_b|testable": "Null B, testable family", "null_b|pr9": "Null B (#9), gate",
         "null_b|pr9_testable": "Null B (#9), testable", "null_b|pr9_by": "Null B (#9), BY", "b2": "**B2 (CAMERA), full**",
         "b2|uncorrected": "B2, no correction (old premise)", "b2|unmatched": "B2, unmatched decoys", "b2|joint": "B2, joint (sides combined)"}
SPATIAL = ["tnbc_wu2021", "allen_brain_atlas", "moffitt2018_hypothalamus"]


def _or(d):
    if not d or d.get("odds_ratio") is None:
        return "–"
    lo, hi = d["ci95"]
    mark = "↑" if lo > 1 else ("↓" if hi < 1 else "")
    return f"{d['odds_ratio']:.2f} [{lo:.2f}, {hi:.2f}]{mark}"


def power_section(lines):
    p = _load("power_a3")
    if not p:
        return
    s, chk = p["summary"], p.get("check_against_power_json", {})
    lines += ["### A3.1 T2 at small effects and B = 5,000", "",
              f"Check: f = 1.5, B = 1,000 reproduces `power.json` (Null A and Null B: calls, power, FDP, one-sided calls) in "
              f"{chk.get('identical', 0)} of {chk.get('checked', 0)} runs.", "",
              "Power / FDR, mean over replicates (FDR 95% t-CI for Null B).", "",
              "| Tool | Setting | Null A (BH) | Null B (BH) | Null B calls (BH) | Null B, tool's rule | B2 (BH) |", "|---|---|---|---|---|---|---|"]
    for t in TOOLS:
        for k, m in s.get(t, {}).items():
            pw = lambda v: f"{m[v]['power'][0]:.2f} / {m[v]['fdr'][0]:.2f}" if v in m else "–"
            lines.append(f"| {t} | {k} | {pw('null_a')} | {pw('null_b')} {_ci(m['null_b']['fdr'][1], 2)} | {m['null_b']['calls']:.1f} | "
                         f"{pw('null_b|tool_rule')} | {pw('b2')} |")
    lines.append("")
    lines += ["**Ablation on T2** (f = 1.5, B = 1,000; power / FDR / one-sided false calls, mean of 5 replicates):", "",
              "| Variant | CellChat | CellPhoneDB |", "|---|---|---|"]
    keys = [k for k in LABEL if any(k in s.get(t, {}).get("f=1.5, B=1000", {}) for t in TOOLS)]
    for v in keys:
        cells = []
        for t in TOOLS:
            m = s.get(t, {}).get("f=1.5, B=1000", {}).get(v)
            cells.append(f"{m['power'][0]:.2f} / {m['fdr'][0]:.2f} / {m['one_sided_false_calls']:.0f}" if m else "–")
        lines.append(f"| {LABEL[v]} | {cells[0]} | {cells[1]} |")
    lines.append("")


def pr9_section(lines):
    rows = []
    for ds in DATASETS:
        e = (_load(f"extend_{ds}") or {}).get("cellphonedb", {}).get("pr9")
        if e:
            v = e["vs_null_b"]
            rows.append(f"| {ds} | {v['calls_b']:,} | {e['calls']['null_b_pr9']:,} | {v['same_call']:.4f} | {v['jaccard']:.3f} | "
                        f"{e['p_identical_in_family']:.3f} | {e['spearman_p_in_family']:.4f} | {e['calls']['null_b_pr9_testable']:,} | "
                        f"{e['calls']['null_b_pr9_by']:,} |")
    if not rows:
        return
    lines += ["### A3.2 CellPhoneDB's Null B as pushed in #9", "",
              "#9's own code: float64 means, `default_rng(0)` permutations; ours: CellPhoneDB's permutation sequence and its float32 means.", "",
              "| Dataset | Calls, ours (BH, gate) | Calls, #9 (BH, gate) | Same call | Jaccard | p identical (family) | Spearman p | #9 testable | #9 BY |",
              "|---|---|---|---|---|---|---|---|---|"] + rows + [""]


def b2_section(lines):
    v = _load("b2_validity")
    lines += ["### A3.3 Null B2: decoy genes with the CAMERA correction", ""]
    if v:
        lines += ["Validity (Tabula Sapiens, 10 replicates). Shuffled labels: side-level type-I error at α = 0.05 (95% t-CI) and "
                  "P(any BH call). Receptor genes shuffled: FDR (95% t-CI). Calibrated = type-I CI reaches ≤ 0.05 and FDR upper bound ≤ 0.05.", "",
                  "| Tissue | Tool | Variant | Type-I (labels shuffled) | P(any call) | FDR (receptor shuffled) | Calibrated |",
                  "|---|---|---|---|---|---|---|"]
        for ds, per in v["results"].items():
            for t, r in per.items():
                for var in ("b2", "b2_uncorrected", "b2_unmatched", "b2_joint"):
                    ls, rs = r["labels_shuffled"][var], r["receptor_shuffled"][var]
                    t1 = ls.get("side_type1", ls.get("row_type1"))
                    cal = r["calibrated"].get(var, {}).get("calibrated")
                    lines.append(f"| {ds} | {t} | {LABEL[var.replace('b2_', 'b2|')]} | {t1[0]:.3f} {_ci(t1[1])} | "
                                 f"{ls['p_any_call']['k']}/{ls['p_any_call']['n']} | {rs['fdr'][0]:.3f} {_ci(rs['fdr'][1])} | "
                                 f"{'–' if cal is None else ('yes' if cal else '**no**')} |")
        lines.append("")
    rows = []
    for ds in DATASETS:
        e = _load(f"extend_{ds}") or {}
        for t in TOOLS:
            b = e.get(t, {}).get("b2")
            if b:
                c = b["calls"]
                rows.append(f"| {ds} | {t} | {b['family']:,} | {b['gated_rows']:,} | {c['p_b2']:,} | {c['p_b2_uncorrected']:,} | "
                            f"{c['p_b2_unmatched']:,} | {c['p_b2_joint']:,} | {b['vs_null_b']['jaccard']:.3f} |")
    if rows:
        lines += ["Real data: BH calls over the tool's family (B2's gate: every subunit in ≥ 10% of its type's cells).", "",
                  "| Dataset | Tool | Family | Passing B2's gate | B2 | No correction | Unmatched | Joint | Jaccard B2 vs Null B |",
                  "|---|---|---|---|---|---|---|---|---|"] + rows + [""]


def ablation_section(lines):
    a = _load("ablation")
    if not a:
        return
    lines += ["### A3.4 Ablation on the five data sets", "",
              "Calls (Jaccard with the full Null B). One part changed at a time.", ""]
    for t in TOOLS:
        dss = [ds for ds in DATASETS if t in a.get(ds, {})]
        lines += [f"**{t}**", "", "| Variant | " + " | ".join(dss) + " |", "|---|" + "---|" * len(dss)]
        for v in LABEL:
            if not any(v in a[ds][t] for ds in dss):
                continue
            cells = [f"{a[ds][t][v]['calls']:,} ({a[ds][t][v]['jaccard_with_null_b']:.2f})" if v in a[ds][t] else "–" for ds in dss]
            lines.append(f"| {LABEL[v]} | " + " | ".join(cells) + " |")
        lines.append("")
    lines += ["Biology of each variant: odds ratio, Null A's calls the variant keeps vs drops (95% bootstrap CI; ↑ above 1, ↓ below 1).", ""]
    for t in TOOLS:
        cols = [(ds, "spatial_kept_vs_dropped", f"{ds} spatial") for ds in SPATIAL] + \
               [(ds, f"D1_{res}", f"{ds} {res}") for ds in DATASETS for res in ("cytosig", "progeny")]
        cols = [c for c in cols if any(c[1] in a.get(c[0], {}).get(t, {}).get(v, {}) for v in LABEL)]
        lines += [f"**{t}**", "", "| Variant | " + " | ".join(c[2] for c in cols) + " |", "|---|" + "---|" * len(cols)]
        for v in LABEL:
            if v == "null_a" or not any(v in a.get(ds, {}).get(t, {}) for ds in DATASETS):
                continue
            lines.append(f"| {LABEL[v]} | " + " | ".join(_or(a.get(ds, {}).get(t, {}).get(v, {}).get(k)) for ds, k, _ in cols) + " |")
        lines.append("")


def amendment3(lines):
    lines += ["## Amendment 3 (an addition; not part of the pre-registered plan)", "",
              f"`AMENDMENT_3.md` (SHA-256 `{amendment_hash(3)[:16]}…`), written after every earlier result and before any of these. "
              "It changes neither verdict above.", ""]
    power_section(lines)
    pr9_section(lines)
    b2_section(lines)
    ablation_section(lines)


STATUS_ROWS = [("A3.1 T2 small effects, B = 5,000", "Amendment 3", ["power_a3"]),
               ("A3.2 #9's CellPhoneDB Null B", "Amendment 3", [f"extend_{d}" for d in DATASETS]),
               ("A3.3 B2 on real data", "Amendment 3", [f"extend_{d}" for d in DATASETS]),
               ("A3.3 B2 validity", "Amendment 3", ["b2_validity"]),
               ("A3.4 ablation", "Amendment 3", ["ablation"])]
