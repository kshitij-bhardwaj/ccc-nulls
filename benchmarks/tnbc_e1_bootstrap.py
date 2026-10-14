"""E1 with uncertainty: AUPRC / odds ratio per tool x null, with 95% bootstrap CIs over the
812 (ligand, target) truth rows, and the paired difference vs the tool's ungated score.

A null that leaves no interaction with p <= 0.05 makes every gated score 0; the AUPRC of an
all-tied ranking is an artefact, so those cells are reported as "no calls".

python benchmarks/tnbc_e1_bootstrap.py full
"""
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import auc, precision_recall_curve

from ccc_nulls import datasets

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "benchmarks" / "results" / "tnbc"
tag = sys.argv[1] if len(sys.argv) > 1 else "full"
NULLS = ["N1", "N2", "N3", "N4", "N5"]
TOOLS = ["CellPhoneDB", "GeometricMean", "CellChat", "Connectome", "NATMI", "SingleCellSignalR", "log2FC"]
NBOOT = 1000

truth = ad.read_h5ad(datasets.path("tnbc_wu2021"), backed="r").uns["ccc_target"][["ligand", "target", "response"]]
truth = truth.astype({"response": int}).reset_index(drop=True)


def joined(rows, score):
    pred = (pd.DataFrame({"ligand": rows.ligand_complex.values, "target": rows.target.values, "score": score})
            .groupby(["ligand", "target"]).score.max().reset_index())
    gt = truth.merge(pred, on=["ligand", "target"], how="left")
    gt["score"] = gt.score.fillna(gt.score.min() - 1e-12 if gt.score.notna().any() else 0)
    return gt.response.to_numpy(), gt.score.to_numpy()


def auprc(y, s):
    p, r, _ = precision_recall_curve(y, s)
    return auc(r, p)


rng = np.random.default_rng(0)
boot = [rng.integers(0, len(truth), len(truth)) for _ in range(NBOOT)]
out = []
for t in TOOLS:
    d = pd.read_csv(RES / f"rows_{t}_{tag}.csv.gz")
    d = d[~d.is_fake]
    y, s0 = joined(d, d.score.to_numpy())
    base = auprc(y, s0)
    base_b = np.array([auprc(y[b], s0[b]) for b in boot])
    out.append(dict(tool=t, null="ungated", auprc=base, lo=np.quantile(base_b, .025), hi=np.quantile(base_b, .975),
                    diff=0.0, diff_lo=0.0, diff_hi=0.0, n_kept=int(len(d)), note=""))
    for nl in NULLS:
        p = d[f"p_{nl}"].to_numpy()
        kept = int((p <= 0.05).sum())
        y1, s1 = joined(d, np.where(p <= 0.05, d.score.to_numpy(), 0.0))
        if kept == 0 or np.unique(s1[s1 > s1.min()]).size == 0:
            out.append(dict(tool=t, null=nl, auprc=np.nan, lo=np.nan, hi=np.nan, diff=np.nan, diff_lo=np.nan,
                            diff_hi=np.nan, n_kept=kept, note="no calls (all scores tied)"))
            continue
        a = auprc(y1, s1)
        ab = np.array([auprc(y1[b], s1[b]) for b in boot])
        db = ab - base_b
        out.append(dict(tool=t, null=nl, auprc=a, lo=np.quantile(ab, .025), hi=np.quantile(ab, .975),
                        diff=a - base, diff_lo=np.quantile(db, .025), diff_hi=np.quantile(db, .975), n_kept=kept,
                        note="CI excludes 0" if (np.quantile(db, .025) > 0 or np.quantile(db, .975) < 0) else ""))
df = pd.DataFrame(out)
df.to_csv(RES / f"e1_bootstrap_{tag}.csv", index=False)
rnd = [auprc(truth.response.to_numpy()[b], np.random.default_rng(i).random(len(truth))) for i, b in enumerate(boot[:200])]
print(f"random baseline AUPRC {np.mean(rnd):.3f} [{np.quantile(rnd, .025):.3f}, {np.quantile(rnd, .975):.3f}]  (positives {truth.response.mean():.3f})")
pd.set_option("display.width", 220)
print(df.assign(auprc=df.apply(lambda r: f"{r.auprc:.3f} [{r.lo:.3f},{r.hi:.3f}]" if pd.notna(r.auprc) else r.note, axis=1),
                diff=df.apply(lambda r: f"{r['diff']:+.3f} [{r.diff_lo:+.3f},{r.diff_hi:+.3f}]" if pd.notna(r['diff']) else "", axis=1))
      [["tool", "null", "auprc", "diff", "n_kept", "note"]].to_string(index=False))
recs = [{k: (None if isinstance(v, float) and v != v else v) for k, v in r.items()} for r in df.to_dict(orient="records")]
json.dump(recs, open(RES / f"e1_bootstrap_{tag}.json", "w"), indent=1, default=float)
