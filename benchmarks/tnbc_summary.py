"""Summarise benchmarks/results/tnbc/results_<tag>.json into tables (markdown + CSV).

python benchmarks/tnbc_summary.py full
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1] / "benchmarks" / "results" / "tnbc"
tag = sys.argv[1] if len(sys.argv) > 1 else "full"
R = json.loads((ROOT / f"results_{tag}.json").read_text())
NULLS = ["N1", "N2", "N3", "N4", "N5"]
NAMES = {"N1": "N1 global perm.", "N2": "N2 patient-strat. perm.", "N3": "N3 matched decoys",
         "N4": "N4 unmatched decoys", "N5": "N5 patient replicates", "N6_native": "N6 LIANA native",
         "ungated": "no null (ungated)"}

rows = []
for t, d in R["tools"].items():
    for key in ["ungated", "N6_native"] + NULLS:
        if key not in d["e1"]:
            continue
        r = dict(tool=t, null=NAMES[key], auprc=d["e1"][key]["auprc"], odds_ratio=d["e1"][key]["odds_ratio"])
        if key in NULLS:
            r.update(auprc_q=d["e1q"][key]["auprc"], or_q=d["e1q"][key]["odds_ratio"],
                     real_calls=d["n_real_calls"][key], real_rate=d["e2"][key]["real_rate"],
                     fake_rate=d["e2"][key]["fake_rate"], fake_p05=d["e2"][key]["fake_p05"],
                     e5_p05=d.get("e5", {}).get(key, {}).get("p05"), e5_calls=d.get("e5", {}).get(key, {}).get("calls"))
        rows.append(r)
df = pd.DataFrame(rows)
df.to_csv(ROOT / f"summary_{tag}.csv", index=False)

pd.set_option("display.width", 200)
print("meta", R["meta"])
print("random baseline", R["random_baseline"])
print("\nvalidation (ours vs LIANA)")
for t, v in R["validation"].items():
    print(f"  {t:18s} max|diff| {v['max_abs_diff']:.1e}  r {v['corr']:.6f}")
print("\nnative vs N1 log-p correlation:", {t: d["native_vs_N1"] for t, d in R["tools"].items() if d["native_vs_N1"] is not None})
print("\nE1 AUPRC (p <= 0.05 gate, Open Problems protocol)")
print(df.pivot(index="tool", columns="null", values="auprc").round(3).to_string())
print("\nE1 AUPRC (q <= 0.05 gate)")
print(df.dropna(subset=["auprc_q"]).pivot(index="tool", columns="null", values="auprc_q").round(3).to_string())
print("\nE2 call rate at q < 0.05, real | fake resource")
e2 = df.dropna(subset=["real_rate"])
print((e2.assign(v=e2.real_rate.map("{:.1%}".format) + " | " + e2.fake_rate.map("{:.1%}".format))
       .pivot(index="tool", columns="null", values="v")).to_string())
print("\nE5 shuffled labels: fraction p < 0.05 (calls at q < 0.05)")
e5 = df.dropna(subset=["e5_p05"])
print((e5.assign(v=e5.e5_p05.map("{:.1%}".format) + " (" + e5.e5_calls.astype(int).astype(str) + ")")
       .pivot(index="tool", columns="null", values="v")).to_string())
