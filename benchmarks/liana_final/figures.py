"""Figures of the benchmark, from the result JSON files only.

    python -m benchmarks.liana_final.figures

Writes PDF to benchmarks/results/liana_final/figures/. Colours: liana6 orange, swap Null B green, 7th/8th Null B blue, the DE filter grey; each
also has its own marker, so no comparison rests on colour alone.
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from benchmarks.liana_nulls.stages import mean_ci

from .pipeline import OUT
from .report2 import DS, cal_summary, power_summary

FIG = OUT / "figures"
INK, MUTE, GRID = "#23373B", "#5E6E72", "#D5DCDD"
STYLE = {"liana6": ("liana6 (baseline)", "#EB811B", "o"), "swap_null_b": ("swap: Null B", "#0B7A2B", "s"),
         "extra_null_b": ("7th/8th: Null B", "#3B6FB6", "D"), "liana6_de": ("filter: both sides DE", "#5E6E72", "^")}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": MUTE, "axes.labelcolor": INK,
                     "xtick.color": MUTE, "ytick.color": MUTE, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "axes.titleweight": "bold", "axes.titlesize": 10, "axes.titlecolor": INK, "legend.frameon": False,
                     "pdf.fonttype": 42})


def J(path):
    return json.loads(path.read_text()) if path.exists() else None


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print("wrote", FIG / f"{name}.pdf")


def calibration():
    """False share in the top 1%, every calibration scenario: F1 (Mammary, Liver), A2.1 (Lung), A2.7 (mouse)."""
    rows = []
    f1 = J(OUT / "calibration.json")
    if f1:
        for t, sides in f1["results"].items():
            for s in ("receptor", "ligand"):
                d = sides[s]
                rows.append((f"{DS[t]}, {s}s", d["false_share_all"][0],
                             {v: tuple(d[v]["false_share_top"]) for v in STYLE if v in d}))
    for part in ("confirm_calibration", "mouse_calibration"):
        j = J(OUT / "amendment2" / f"{part}.json")
        if j:
            S = cal_summary(j["runs"])
            for (d, setting, s), m in sorted(S.items()):
                if setting == "default":
                    rows.append((f"{DS[d]}, {s}s" + (" (held out)" if d == "tabula_sapiens_lung" else ""), m["all"][0],
                                 {v: (m[v][0], m[v][1]) for v in STYLE if v in m}))
    if not rows:
        return
    fig, ax = plt.subplots(figsize=(6.4, 0.42 * len(rows) + 1.0))
    y = np.arange(len(rows))[::-1]
    off = dict(zip(STYLE, np.linspace(0.27, -0.27, len(STYLE))))
    for yi, (lab, allv, vals) in zip(y, rows):
        ax.plot([allv], [yi], marker="|", ms=14, color=MUTE, mew=1.5)
        for v, (m, ci) in vals.items():
            name, col, mk = STYLE[v]
            ax.errorbar(m, yi + off[v], xerr=[[m - ci[0]], [ci[1] - m]], fmt=mk, color=col, ms=5, elinewidth=1, capsize=2)
    ax.set_yticks(y, [r[0] for r in rows])
    ax.set_xlabel("share of the top 1% whose shuffled side is null  (lower is better)")
    ax.set_title("Calibration: one-sided rows among the top-ranked", loc="left")
    handles = [plt.Line2D([], [], marker=STYLE[v][2], color=STYLE[v][1], ls="", label=STYLE[v][0]) for v in STYLE]
    handles.append(plt.Line2D([], [], marker="|", ms=10, mew=1.5, color=MUTE, ls="", label="share among all rows"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12 - 0.3 / len(rows)), ncol=3, fontsize=8)
    save(fig, "calibration")


def power():
    """A2.2: AUROC planted vs one-sided and power in the top 1%, by effect size and class."""
    j = J(OUT / "amendment2" / "hard_power.json")
    if not j:
        return
    P = power_summary(j["runs"])
    keys = sorted(P, key=lambda k: k[3])
    fs = [k[3] for k in keys]
    classes = [("single", "single pairs (30)"), ("complex", "complexes (10)"), ("unequal", "unequal sides (10)")]
    fig, axes = plt.subplots(2, 3, figsize=(9.6, 5.2), sharex=True, sharey="row")
    for c, (cls, title) in enumerate(classes):
        for r, (m, lab) in enumerate((("auroc", "AUROC, planted vs one-sided"), ("power_top", "power in the top 1%"))):
            ax = axes[r, c]
            for v, (name, col, mk) in STYLE.items():
                if v not in P[keys[0]]:
                    continue
                mu = [P[k][v][cls][m][0] for k in keys]
                lo = [P[k][v][cls][m][1][0] for k in keys]
                hi = [P[k][v][cls][m][1][1] for k in keys]
                ax.plot(fs, mu, marker=mk, color=col, ms=5, lw=1.6, label=name)
                ax.fill_between(fs, lo, hi, color=col, alpha=0.12, lw=0)
            if r == 0:
                ax.set_title(title, loc="left")
            if c == 0:
                ax.set_ylabel(lab)
            if r == 1:
                ax.set_xlabel("fold change f")
                ax.set_xticks(fs)
    axes[0, 0].legend(loc="lower right", fontsize=7.5)
    fig.suptitle("Power on planted edges (TS Mammary background)", x=0.01, ha="left", fontweight="bold", color=INK)
    fig.tight_layout()
    save(fig, "power")


def pvalues():
    """A2.6: Q-Q of the tool-level p on false rows (shuffled side), pooled over the calibration reruns."""
    j = J(OUT / "amendment2" / "recalibration.json")
    if not j:
        return
    runs = [r for r in j["runs"].values() if "p_calibration" in r]
    cols = {"null_a": ("Null A (LIANA+'s p)", "#EB811B", "-"), "null_b": ("Null B", "#0B7A2B", "-"), "b2": ("B2", "#3B6FB6", "--")}
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6), sharey=True)
    u = np.linspace(0, 1, 101)
    for ax, tool in zip(axes, ("CellPhoneDB", "CellChat")):
        ax.plot([0, 1], [0, 1], color=MUTE, lw=1, ls=":")
        for n, (lab, col, ls) in cols.items():
            q = np.mean([r["p_calibration"][tool][n]["quantiles"] for r in runs if r["p_calibration"][tool][n]["quantiles"]], axis=0)
            ax.plot(u, q, color=col, ls=ls, lw=1.8, label=lab)
        ax.set_title(tool, loc="left")
        ax.set_xlabel("uniform quantile")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    axes[0].set_ylabel("p-value quantile on false rows")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle("Below the diagonal: anti-conservative on one-sided rows", x=0.01, ha="left", fontweight="bold", color=INK)
    fig.tight_layout()
    save(fig, "pvalues")


def donor():
    """A2.3 (a): share of p <= 0.05 with labels shuffled within donor."""
    j = J(OUT / "amendment2" / "donor_null.json")
    if not j:
        return
    runs = list(j["runs"].values())
    ds = list(dict.fromkeys(r["dataset"] for r in runs))
    cats = [(d, t) for d in ds for t in ("CellPhoneDB", "CellChat")]
    series = [("null_a", "global", "Null A, global", "#EB811B", ""), ("null_a", "donor", "Null A, within donor", "#EB811B", "///"),
              ("null_b", "global", "Null B, global", "#0B7A2B", ""), ("null_b", "donor", "Null B, within donor", "#0B7A2B", "///")]
    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(cats)), 3.6))
    x = np.arange(len(cats))
    w = 0.8 / len(series)
    for i, (n, sch, lab, col, h) in enumerate(series):
        vals = [mean_ci([r[t][f"{n}_le_0.05"] for r in runs if r["dataset"] == d and r["scheme"] == sch] or [np.nan])[0] for d, t in cats]
        ax.bar(x - 0.4 + w * (i + 0.5), vals, w * 0.92, color=col if not h else "white", edgecolor=col, hatch=h, lw=1, label=lab)
    ax.axhline(0.05, color=INK, ls="--", lw=1)
    ax.text(len(cats) - 0.5, 0.05, "nominal 0.05", ha="right", va="bottom", fontsize=8, color=INK)
    ax.set_xticks(x, [f"{DS[d]}\n{t}" for d, t in cats], fontsize=8)
    ax.set_ylabel("share of rows with p ≤ 0.05")
    ax.set_title("Labels shuffled within donor: no cell-type signal is left", loc="left")
    ax.legend(ncol=2, fontsize=8, loc="upper left")
    save(fig, "donor")


def biology():
    """A2.8: odds ratios of every F3 test, with the detectable range."""
    j = J(OUT / "amendment2" / "biology.json")
    if not j or "tests" not in j:
        return
    T = j["tests"]
    variants = [("swap_null_b", "swap: Null B"), ("extra_null_b", "7th/8th: Null B")]
    tests = list(dict.fromkeys(k.split(" | ")[1] for k in T if k.startswith("swap_null_b |")))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 0.33 * len(tests) + 1.2), sharey=True)
    y = np.arange(len(tests))[::-1]
    for ax, (v, title) in zip(axes, variants):
        for yi, t in zip(y, tests):
            r = T.get(f"{v} | {t}")
            if not r or r["odds_ratio"] is None:
                continue
            lo, hi = r["detectable"]
            ax.fill_betweenx([yi - 0.35, yi + 0.35], lo, hi, color=GRID, lw=0)
            col = "#0B7A2B" if r.get("q_all", 1) < 0.05 and r["odds_ratio"] > 1 else ("#B33A3A" if r.get("q_all", 1) < 0.05 else INK)
            ax.errorbar(r["odds_ratio"], yi, xerr=[[r["odds_ratio"] - r["ci95"][0]], [r["ci95"][1] - r["odds_ratio"]]],
                        fmt="o", color=col, ms=4, elinewidth=1, capsize=2)
        ax.axvline(1, color=MUTE, lw=1)
        ax.set_xscale("log")
        ax.set_title(title, loc="left")
        ax.set_xlabel("odds ratio, kept vs dropped (log scale)")
    axes[0].set_yticks(y, [t.replace("_", " ") for t in tests], fontsize=8)
    fig.text(0.01, -0.02, "Grey band: odds ratios the test could not detect with 80% power. Coloured: BH q < 0.05 across all tests "
             "(green above 1, red below 1).", fontsize=7.5, color=MUTE, ha="left", va="top")
    fig.tight_layout()
    save(fig, "biology")


def main():
    for f in (calibration, power, pvalues, donor, biology):
        f()


if __name__ == "__main__":
    main()
