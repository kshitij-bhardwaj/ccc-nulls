"""benchmarks/results/liana_final/amendment6/README.md from the AMENDMENT_6 result files."""
from __future__ import annotations

import json

import numpy as np

from .amendment6 import NULLS4, OUT6, SHOWN, SWAP2, SWAP3
from .pipeline import BASELINE, OUT
from .report3 import mci, table
from .report4 import crit1, crit2, word

EARLIER = {"calibration": OUT / "amendment3/calibration.json", "confirm_calibration": OUT / "amendment4/confirm_calibration.json",
           "mouse_calibration": OUT / "amendment4/mouse_calibration.json", "power": OUT / "amendment3/power.json",
           "hard_power": OUT / "amendment3/hard_power.json", "confirm_power": OUT / "amendment4/confirm_power.json"}
SETTINGS = [("calibration", "tabula_sapiens_mammary", "Mammary"), ("calibration", "tabula_sapiens_liver", "Liver"),
            ("confirm_calibration", "tabula_sapiens_lung", "Lung (held out)"), ("mouse_calibration", "allen_brain_atlas", "Allen brain")]


def J(path):
    return json.loads(path.read_text()) if path.exists() else None


def check(part):
    """Runs in both files whose liana6 and earlier variants are identical (calibration: false share; power: AUROC)."""
    new, old = J(OUT6 / f"{part}.json"), J(EARLIER[part])
    if not (new and old):
        return None
    both = [t for t in new["runs"] if t in old["runs"]]
    if "power" in part:
        get = lambda r, v: (r[v]["single"]["power_top"], r[v]["single"]["auroc"])
    else:
        get = lambda r, v: r["false_share_top"][v]
    n = sum(all(get(new["runs"][t], v) == get(old["runs"][t], v) for v in [BASELINE, *SWAP2]) for t in both)
    return f"{n} of {len(both)}"


def write():
    f = lambda n: J(OUT6 / f"{n}.json")
    L = ["# AMENDMENT_6 results: Connectome's p-value replaced too", "",
         "Plan: [`AMENDMENT_6.md`](../../../liana_final/AMENDMENT_6.md). `swap_<n>`: CellPhoneDB's and CellChat's p "
         "replaced by null n's; `swap3_<n>`: Connectome's specificity too. Same runs, so the columns are paired. "
         "Intervals: 95% t-intervals over replicates.", ""]

    L += ["## Calibration: false share of the top 1% (lower is better)", ""]
    rows = []
    for part, name, short in SETTINGS:
        cal = f(part)
        for s in ("receptor", "ligand"):
            rs = [r for r in (cal or {"runs": {}})["runs"].values() if r["dataset"] == name and r["scenario"] == s
                  and r.get("setting", "default") == "default"]
            if rs:
                rows.append([f"{short}, {s}s", len(rs)] + [mci([r["false_share_top"][v] for r in rs]) for v in SHOWN])
    L += table(["Setting", "reps", *SHOWN], rows)

    cc, cp = f("confirm_calibration"), f("confirm_power")
    L += ["## Held-out Lung: AMENDMENT_1's rule", ""]
    rows = []
    for v in [*SWAP2, *SWAP3]:
        c1 = crit1(cc, "tabula_sapiens_lung", v) if cc else None
        c2 = crit2(cp, v) if cp else None
        out = None if None in (c1, c2) else "confirmed" if c1 and c2 else \
            "not confirmed (" + ", ".join(n for n, c in (("calibration", c1), ("power", c2)) if not c) + ")"
        rows.append([v, word(c1, "holds", "fails"), word(c2, "holds", "fails"), out or "pending"])
    L += table(["Variant", "1 calibration", "2 power", "Outcome"], rows)

    for part, title in (("power", "Power, F2 design"), ("hard_power", "Power, small effects (A2.2)"), ("confirm_power", "Power, held-out Lung")):
        pw = f(part)
        if not pw:
            continue
        L += [f"## {title}: power in the top 1% / AUROC, single pairs", ""]
        rows = []
        for e in sorted({r["effect"] for r in pw["runs"].values()}):
            rs = [r for r in pw["runs"].values() if r["effect"] == e]
            rows.append([f"f = {e}", len(rs)] + [f"{np.mean([r[v]['single']['power_top'] for r in rs]):.2f} / "
                                                f"{np.nanmean([r[v]['single']['auroc'] for r in rs]):.3f}" for v in SHOWN])
        L += table(["Effect", "reps", *SHOWN], rows)

    cal = f("calibration")
    if cal:
        L += ["## Connectome's own p on false rows (Mammary and Liver; share ≤ 0.05 / BH discoveries, false share)", ""]
        rs = [r for r in cal["runs"].values() if "Connectome" in r.get("p_calibration", {})]
        rows = []
        for n in NULLS4:
            x = [r["p_calibration"]["Connectome"][n] for r in rs]
            fs = [y["bh_false_share"] for y in x if y["bh_false_share"] is not None]
            rows.append([n, len(rs), mci([y["at_most"]["0.05"] for y in x]), f"{np.mean([y['bh_discoveries'] for y in x]):,.0f}",
                         f"{np.mean(fs):.3f}" if fs else "–"])
        L += table(["Null", "runs", "p ≤ 0.05", "BH discoveries", "false share of discoveries"], rows)

    real = f("real")
    if real:
        L += ["## Real labels: calls (specificity_rank ≤ 0.05) and Jaccard with liana6", ""]
        L += table(["Data set", *SHOWN], [[n] + [f"{r['calls'][v]['calls']:,} ({r['calls'][v]['jaccard_with_liana6']:.2f})" for v in SHOWN]
                                          for n, r in real["runs"].items()])
    bio = f("biology")
    if bio:
        t = bio["tests"]
        L += ["## Biology (A2.8): spatial odds ratio, kept vs dropped [95% CI], q across all", ""]
        rows = []
        for v in ["swap_null_b", "swap_b2", "swap_t2", *SWAP3]:
            rows.append([v] + [(lambda r: f"{r['odds_ratio']:.2f} [{r['ci95'][0]:.2f}, {r['ci95'][1]:.2f}] q {r['q_all']:.2f}"
                                if r.get("odds_ratio") is not None else "–")(t.get(f"{v} | {d} spatial", {}))
                               for d in ("tnbc_wu2021", "allen_brain_atlas", "moffitt2018_hypothalamus")])
        L += table(["Variant", "TNBC", "Allen brain", "hypothalamus"], rows)
        sig = [k for k, r in t.items() if r.get("q_all") is not None and r["q_all"] < 0.05]
        L += [f"Tests with q (all) < 0.05: {', '.join(sig) if sig else 'none'}.", ""]

    L += ["## Checks: liana6 and the earlier variants reproduce AMENDMENTS 3 and 4", ""]
    L += table(["File", "runs identical"], [[p, check(p) or "pending"] for p in EARLIER])
    OUT6.mkdir(parents=True, exist_ok=True)
    (OUT6 / "README.md").write_text("\n".join(L) + "\n")
    print("wrote", OUT6 / "README.md")
