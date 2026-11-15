"""The paper's three figures. Each chart puts all four tools (CellChat, CellPhoneDB, Connectome, NATMI) and all four
tests (JLP = the tools' own test, BPT, MDT, LOO-t) on one set of axes, instead of one panel per tool.

    python figures/make_figures.py            (from code/; writes figures/output/fig{1,2,3}_*.pdf and .png)

Reads the committed result files (benchmarks/results/liana_nulls and liana_final) and figures/final_results.json
(the plan's Benchmark-1 numbers, copied from the PLAN report). Nothing is recomputed, so it runs in seconds. Benchmark 2 is read from AMENDMENT_8's
files, which recompute every earlier variant in the same runs (checked: identical to AMENDMENTS 3, 4 and 6).
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
OUT = HERE / "output"
REPO = HERE.parent                                  # code/
R1, R2 = REPO / "benchmarks/results/liana_nulls", REPO / "benchmarks/results/liana_final"
R4, R5, R7, R8 = (R2 / f"amendment{i}" for i in (4, 5, 7, 8))
PLAN = json.loads((HERE / "final_results.json").read_text())

COL = {"JLP": "#EB811B", "BPT": "#0B7A2B", "MDT": "#8C6D31", "LOO-t": "#7B3FA0"}
LIGHT = {"JLP": "#F5C08A", "BPT": "#8FCF9F"}
TOOLS = ["CellChat", "CellPhoneDB", "Connectome", "NATMI"]
MARK = {"CellChat": "o", "CellPhoneDB": "s", "Connectome": "^", "NATMI": "D"}
# Benchmark 2: the steps along the x axis, i.e. which tools' p-values the test replaces in LIANA+'s consensus
STEPS = ["as shipped", "CellChat +\nCellPhoneDB", "+ Connectome", "+ NATMI"]
VAR = {"BPT": ("null_b", ["swap_null_b", "swap3_null_b", "swap4_null_b"]),
       "MDT": ("b2", ["swap_b2", "swap3_b2", "swap4_b2"]),
       "LOO-t": ("t2", ["swap_t2", "swap3_t2", "swap4_t2"])}

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#5E6E72",
                     "axes.labelcolor": "#23373B", "xtick.color": "#5E6E72", "ytick.color": "#5E6E72",
                     "axes.titlesize": 8.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
                     "font.family": "DejaVu Sans", "legend.frameon": False, "pdf.fonttype": 42})


def J(p):
    return json.loads(Path(p).read_text())


def ci(x):
    x = np.asarray(x, float)
    h = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0
    return x.mean(), h


def style(ax, ylabel=None):
    ax.grid(axis="y", color="#E3E7E8", lw=0.5); ax.set_axisbelow(True)
    if ylabel:
        ax.set_ylabel(ylabel)


def tag(ax, letter):
    ax.text(-0.13, 1.06, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def grouped(ax, cats, series, width=0.82):
    """Grouped bars. series: [(name, colour, means, (lo, hi) whiskers or None)]; NaN bars are skipped."""
    x, w = np.arange(len(cats)), width / len(series)
    for j, (name, col, v, e) in enumerate(series):
        xs = x + (j - (len(series) - 1) / 2) * w
        v = np.asarray(v, float)
        ax.bar(xs, v, w, color=col, label=name, edgecolor="white", linewidth=0.3)
        if e is not None:
            lo, hi = np.asarray(e[0], float), np.asarray(e[1], float)
            ax.errorbar(xs, v, yerr=[v - lo, hi - v], fmt="none", lw=0.6, capsize=1.2, ecolor="#23373B")
    ax.set_xticks(x); ax.set_xticklabels(cats)


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight"); fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=220)
    plt.close(fig); print("wrote", name)


# ------------------------------------------------------------------ Benchmark 1: each tool's own p-values
def b1_known_nulls():
    """FDR at BH 5% when one side's genes are shuffled (every row on that side is false), per tool and test,
    over 2 tissues x 2 sides. CellChat, CellPhoneDB: the plan's runs (MDT, LOO-t: AMENDMENT_4); Connectome, NATMI:
    inside LIANA+ (AMENDMENTS 6 and 8). No discovery counts as FDR 0."""
    a4 = J(R1 / "amendment4_calibration.json")["results"]
    out = {t: {n: [] for n in COL} for t in TOOLS}
    for side, key in (("receptor", "calibration_receptor_fdr"), ("ligand", "calibration_ligand_fdr")):
        plan = PLAN[key]
        for i, (t, s) in enumerate([(t, s) for t in ("cellchat", "cellphonedb")
                                     for s in ("tabula_sapiens_mammary", "tabula_sapiens_liver")]):
            name = {"cellchat": "CellChat", "cellphonedb": "CellPhoneDB"}[t]
            out[name]["JLP"].append(plan["null_a"][i]); out[name]["BPT"].append(plan["null_b"][i])
            out[name]["MDT"].append(a4[s][t]["summary"]["mdt"][side]["fdr"][0])
            out[name]["LOO-t"].append(a4[s][t]["summary"]["loo_t"][side]["fdr"][0])
    for tool, f in (("Connectome", R2 / "amendment6/calibration.json"), ("NATMI", R8 / "calibration.json")):
        rr = list(J(f)["runs"].values())
        for side in ("receptor", "ligand"):
            for d in ("tabula_sapiens_mammary", "tabula_sapiens_liver"):
                for k, n in (("null_a", "JLP"), ("null_b", "BPT"), ("b2", "MDT"), ("t2", "LOO-t")):
                    out[tool][n].append(float(np.mean([r["p_calibration"][tool][k]["bh_false_share"] or 0.0 for r in rr
                                                       if r["dataset"] == d and r["scenario"] == side])))
    return out


def b1_power():
    """(power, FDR) on 30 planted pairs, f = 1.5, BH 5%: CellChat and CellPhoneDB from the plan (LOO-t: AMENDMENT_4),
    Connectome (AMENDMENT_7) and NATMI (AMENDMENT_8) inside LIANA+, JLP there with (1 + #)/(B + 1)."""
    p, a4 = PLAN["power_f15"], J(R1 / "amendment4_power.json")["summary"]
    out = {}
    for i, t in enumerate(("CellChat", "CellPhoneDB")):
        out[t] = {"JLP": (p["null_a"]["power"][i], p["null_a"]["fdr"][i]),
                  "BPT, B = 1,000": (p["null_b_B1000"]["power"][i], p["null_b_B1000"]["fdr"][i]),
                  "BPT": (p["null_b_B5000"]["power"][i], p["null_b_B5000"]["fdr"][i]),
                  "MDT": (p["null_b2"]["power"][i], p["null_b2"]["fdr"][i]),
                  "LOO-t": (a4[t.lower()]["power"][0], a4[t.lower()]["fdr"][0])}
    for t, f in (("Connectome", R7 / "tool_power.json"), ("NATMI", R8 / "tool_power.json")):
        rr = list(J(f)["runs"].values())
        m = lambda n, B, k: float(np.mean([r[t][n][k] for r in rr if r["n_perms"] == B]))
        out[t] = {lab: (m(n, B, "power"), m(n, B, "fdp")) for lab, n, B in
                  (("JLP", "null_a_plus1", 1000), ("BPT, B = 1,000", "null_b", 1000), ("BPT", "null_b", 5000),
                   ("MDT", "b2", 1000), ("LOO-t", "t2", 1000))}
    return out


def fig1():
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.6, 2.45), gridspec_kw={"width_ratios": [1.35, 1]})
    kn = b1_known_nulls()
    series = [(n, c, [100 * np.mean(kn[t][n]) for t in TOOLS],
               ([100 * min(kn[t][n]) for t in TOOLS], [100 * max(kn[t][n]) for t in TOOLS])) for n, c in COL.items()]
    grouped(a, TOOLS, series)
    a.axhline(5, ls="--", lw=0.7, color="#23373B"); a.text(3.45, 7, "5%", fontsize=6.5, color="#23373B", ha="right")
    for i, t in enumerate(TOOLS):            # bars too low to see: print their value
        for j, (n, c, v, _) in enumerate(series):
            if v[i] < 0.5:
                a.text(i + (j - 1.5) * 0.205, 1.2, "0", ha="center", fontsize=5.5, color=c)
    a.set_ylim(0, 105); style(a, "FDR at BH 5% (%)"); a.tick_params(axis="x", labelsize=7)
    a.set_title("Known nulls: one side's genes shuffled", pad=14)
    a.legend(ncol=4, fontsize=7, loc="lower left", bbox_to_anchor=(0, 0.98), handlelength=1, columnspacing=1)
    tag(a, "a")

    pw = b1_power()
    jit = {"CellChat": (-0.012, 0.0), "CellPhoneDB": (0.012, 0.0), "Connectome": (-0.012, -0.02), "NATMI": (0.012, 0.02)}
    for t in TOOLS:
        for lab, (pv, fd) in pw[t].items():
            n = lab.split(",")[0]
            dx, dy = jit[t]
            b.scatter(100 * (fd + dx), 100 * (pv + dy), marker=MARK[t], s=26, color=COL[n] if "1,000" not in lab else "white",
                      edgecolor=COL[n], linewidth=0.9, zorder=3)
    b.set_xlim(-6, 104); b.set_ylim(-6, 106)
    b.set_xlabel("FDR among calls (%)"); style(b, "power (%)")
    b.text(2, 101, "ideal", fontsize=6.5, color="#5E6E72", va="bottom")
    b.set_title("30 planted pairs, fold change 1.5", pad=14)
    h = [Line2D([], [], ls="", marker=MARK[t], color="#5E6E72", ms=4.5, label=t) for t in TOOLS] + \
        [Line2D([], [], ls="", marker="o", mfc="white", mec=COL["BPT"], ms=4.5, label="BPT, 1,000 perm.")]
    b.legend(handles=h, fontsize=6.3, loc="center right", handletextpad=0.2, borderaxespad=0.2)
    tag(b, "b")
    fig.tight_layout(w_pad=2.2)
    save(fig, "fig1_tool_level")


# ------------------------------------------------------------------ Benchmark 2: LIANA+'s consensus
def b2_runs(f, d):
    return [r for r in J(R8 / f)["runs"].values() if r.get("dataset") == d and r.get("setting", "default") == "default"]


def steps(rr, metric, test):
    """Mean and 95% CI of `metric` over runs rr at each step: baseline, then the test replacing 2, 3, 4 tools' p."""
    keys = ["liana6"] + VAR[test][1]
    return [ci([metric(r, k) for r in rr]) for k in keys]


def fig2():
    fig = plt.figure(figsize=(6.6, 4.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0.75, wspace=0.42)
    gb = fig.add_gridspec(2, 2, height_ratios=[1, 1], hspace=0.75, wspace=0.36)
    panels = [("Mammary", "calibration.json", "tabula_sapiens_mammary"), ("Liver", "calibration.json", "tabula_sapiens_liver"),
              ("Lung (held out)", "confirm_calibration.json", "tabula_sapiens_lung"),
              ("Mouse brain", "mouse_calibration.json", "allen_brain_atlas")]
    x = np.arange(4)
    for i, (title, f, d) in enumerate(panels):
        ax = fig.add_subplot(gs[0, i])
        rr = b2_runs(f, d)           # 3 replicates x receptor/ligand shuffled
        for n in VAR:
            v = steps(rr, lambda r, k: 100 * r["false_share_top"][k], n)
            ax.errorbar(x, [m for m, _ in v], yerr=[h for _, h in v], color=COL[n], marker="o", ms=3.2, lw=1.4,
                        capsize=1.5, label=n)
        ax.plot(0, np.mean([100 * r["false_share_top"]["liana6"] for r in rr]), "o", color=COL["JLP"], ms=5, zorder=4)
        ax.set_ylim(30, 100); ax.set_xticks(x); ax.set_xticklabels(["0", "2", "3", "4"])
        ax.set_title(title, fontsize=8); style(ax, "false share of top 1% (%)" if i == 0 else None)
        if i == 0:
            tag(ax, "a")
    fig.text(0.5, 0.505, "number of tools whose p-value the test replaces (0 = LIANA+ as shipped, JLP; "
             "2 = CellChat, CellPhoneDB; 3 = + Connectome; 4 = + NATMI)", ha="center", fontsize=7, color="#23373B")

    hp = [r for k, r in J(R8 / "hard_power.json")["runs"].items() if "|f1.25|" in k]
    for j, (metric, lab, lim) in enumerate([("power_top", "power in top 1% (%)", (50, 102)),
                                            ("auroc", "AUROC", (0.78, 1.0))]):
        ax = fig.add_subplot(gb[1, j])
        scale = 100 if metric == "power_top" else 1
        for n in VAR:
            v = steps(hp, lambda r, k: scale * r[k]["single"][metric], n)
            ax.errorbar(x, [m for m, _ in v], yerr=[h for _, h in v], color=COL[n], marker="o", ms=3.2, lw=1.4,
                        capsize=1.5, label=n)
        ax.plot(0, np.mean([scale * r["liana6"]["single"][metric] for r in hp]), "o", color=COL["JLP"], ms=5, zorder=4,
                label="JLP (LIANA+ as shipped)")
        ax.set_xticks(x); ax.set_xticklabels(["as shipped", "2 tools", "+ Connectome", "+ NATMI"], fontsize=7)
        ax.set_ylim(*lim); style(ax, lab)
        ax.set_title("Planted pairs: found in the top 1%" if j == 0 else "Planted pairs vs one-sided rows", fontsize=8)
        if j == 0:
            tag(ax, "b")
    h, l = ax.get_legend_handles_labels()
    order = [l.index("JLP (LIANA+ as shipped)")] + [l.index(n) for n in VAR]
    fig.legend([h[i] for i in order], [l[i] for i in order], ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=7.5)
    save(fig, "fig2_consensus")


# ------------------------------------------------------------------ p-value validity, donors, robustness, biology
def fig3():
    fig, axs = plt.subplots(2, 2, figsize=(6.6, 4.9), gridspec_kw={"hspace": 0.62, "wspace": 0.26})
    # (a) one side shuffled inside LIANA+: share of false rows with p <= 0.05 (Mammary and Liver, 12 runs)
    ax = axs[0, 0]
    rr = list(J(R8 / "calibration.json")["runs"].values())
    series = []
    for k, n in (("null_a", "JLP"), ("null_b", "BPT"), ("b2", "MDT"), ("t2", "LOO-t")):
        v = [ci([100 * r["p_calibration"][t][k]["at_most"]["0.05"] for r in rr]) for t in TOOLS]
        series.append((n, COL[n], [m for m, _ in v], ([m - h for m, h in v], [m + h for m, h in v])))
    grouped(ax, TOOLS, series)
    ax.axhline(5, ls="--", lw=0.7, color="#23373B"); ax.set_ylim(0, 50); ax.tick_params(axis="x", labelsize=5.6)
    style(ax, "false rows with p ≤ 0.05 (%)"); ax.set_title("One side shuffled, inside LIANA+")
    ax.legend(ncol=4, fontsize=6.5, loc="upper right", handlelength=1, columnspacing=0.8)
    tag(ax, "a")

    # (b) labels shuffled within donor (no cell-type signal): rows with p <= 0.05, mean over 4 data sets (min-max)
    ax = axs[0, 1]
    a2, a4, a5 = J(R2 / "amendment2/donor_null.json")["runs"], J(R4 / "donor_null.json")["runs"], J(R5 / "donor_b2.json")["runs"]
    a7, a8 = J(R7 / "donor_null.json")["runs"], J(R8 / "donor_null.json")["runs"]
    ds = ("tabula_sapiens_lung", "tnbc_wu2021", "tabula_sapiens_mammary", "tabula_sapiens_liver")

    def donor(t, src, scheme, key):
        s_ = {"Connectome": a7, "NATMI": a8}.get(t, src)
        sc = "global" if t in ("Connectome", "NATMI") and scheme is None else scheme
        return [100 * float(np.mean([r[t][key] for r in s_.values() if r["dataset"] == d and (sc is None or r.get("scheme") == sc)]))
                for d in ds]
    rows = [("JLP, global", COL["JLP"], a2, "global", "null_a_le_0.05"), ("JLP, within donor", LIGHT["JLP"], a2, "donor", "null_a_le_0.05"),
            ("BPT, global", COL["BPT"], a2, "global", "null_b_le_0.05"), ("BPT, within donor", LIGHT["BPT"], a2, "donor", "null_b_le_0.05"),
            ("MDT", COL["MDT"], a5, None, "b2_le_0.05"), ("LOO-t", COL["LOO-t"], a4, "global", "t2_le_0.05")]
    series = []
    for name, c, src, sc, key in rows:
        v = [donor(t, src, sc, key) for t in TOOLS]
        series.append((name, c, [np.mean(x) for x in v], ([min(x) for x in v], [max(x) for x in v])))
    grouped(ax, TOOLS, series)
    ax.axhline(5, ls="--", lw=0.7, color="#23373B"); ax.set_ylim(0, 75); ax.tick_params(axis="x", labelsize=5.6)
    style(ax, "rows with p ≤ 0.05 (%)"); ax.set_title("Labels shuffled within donor")
    ax.legend(ncol=3, fontsize=6, loc="upper right", handlelength=1, columnspacing=0.8)
    tag(ax, "b")

    # (c) sensitivity: false share of the top 1% minus the baseline's, per run (settings, permutations, mouse)
    ax = axs[1, 0]
    s4 = [r for f in ("sensitivity_calibration.json", "sensitivity_permutations_calibration.json", "mouse_calibration.json")
          for r in J(R4 / f)["runs"].values()]
    s7 = [r for f in ("sensitivity_calibration.json", "sensitivity_permutations_calibration.json") for r in J(R7 / f)["runs"].values()]
    s8 = [r for f in ("sensitivity_calibration.json", "sensitivity_permutations_calibration.json") for r in J(R8 / f)["runs"].values()]
    rng = np.random.default_rng(0)
    for i, n in enumerate(VAR):
        for j, (src, k) in enumerate(zip((s4, s7, s8), VAR[n][1])):
            dv = [100 * (r["false_share_top"][k] - r["false_share_top"]["liana6"]) for r in src]
            xs = i + (j - 1) * 0.26 + rng.uniform(-0.07, 0.07, len(dv))
            ax.scatter(xs, dv, s=7, color=COL[n], alpha=[1.0, 0.65, 0.4][j], linewidth=0)
    ax.axhline(0, color="#23373B", lw=0.8)
    ax.set_xticks([i + d for i in range(3) for d in (-0.26, 0, 0.26)])
    ax.set_xticklabels([x for n in VAR for x in ("2", f"3\n{n}", "4")])     # tools replaced, under each test
    ax.tick_params(axis="x", labelsize=6.5)
    style(ax, "top-1% false share minus\nbaseline (points)"); ax.set_title("Sensitivity runs (below 0: better)")
    ax.set_xlabel("tools whose p-value the test replaces", fontsize=6.5)
    tag(ax, "c")

    # (d) spatial co-localisation among the baseline's calls, kept vs dropped (odds ratio, 95% CI; grey: not
    # detectable with 80% power)
    ax = axs[1, 1]
    bio = J(R8 / "biology.json")["tests"]
    dsb = [("tnbc_wu2021", "TNBC\n(Visium)"), ("allen_brain_atlas", "Allen brain"), ("moffitt2018_hypothalamus", "Hypothalamus\n(MERFISH)")]
    for i, (d, _) in enumerate(dsb):
        det = np.array([bio[f"{VAR[n][1][0]} | {d} spatial"]["detectable"] for n in VAR])
        ax.fill_between([i - 0.45, i + 0.45], np.log10(det[:, 0].max()), np.log10(det[:, 1].min()), color="#E3E7E8", lw=0)
        for a, n in enumerate(VAR):
            for j, k in enumerate(VAR[n][1]):
                r = bio[f"{k} | {d} spatial"]
                xx = i + (a - 1) * 0.28 + (j - 1) * 0.08
                ax.plot([xx, xx], np.log10(r["ci95"]), color=COL[n], lw=0.9, alpha=[1, 0.7, 0.45][j])
                ax.plot(xx, np.log10(r["odds_ratio"]), "o", ms=2.6, color=COL[n], alpha=[1, 0.7, 0.45][j])
    ax.axhline(0, color="#23373B", lw=0.8)
    ax.set_xticks(range(3)); ax.set_xticklabels([n for _, n in dsb], fontsize=6.8)
    ax.set_yticks(np.log10([0.25, 0.5, 1, 2, 4])); ax.set_yticklabels(["0.25", "0.5", "1", "2", "4"])
    style(ax, "odds ratio, kept vs dropped"); ax.set_title("Spatial co-localisation of kept calls")
    tag(ax, "d")
    save(fig, "fig3_robustness")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    only = set(sys.argv[1:])
    for f in (fig1, fig2, fig3):
        if not only or f.__name__ in only:
            f()
