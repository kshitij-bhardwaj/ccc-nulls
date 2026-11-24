"""NicheNet footprints for NATMI, CellChat and CellPhoneDB (NOT pre-registered; an addition after AMENDMENT_5).

    python -m benchmarks.liana_nulls.nichenet_eval [--tool natmi cellchat cellphonedb] [--dataset ...]

Step 1  truth:   NicheNet as a third downstream truth next to CytoSig and PROGENy (Amendment 1, A2): among Null A's
                 calls with s != t, odds that the receiver shows the ligand's footprint, calls Null B keeps vs drops
                 (D1; D2 by limiting partner). Same scoring and bootstrap as the other two resources.
Step 2  p-value: footprint p(L, t) (nichenet.footprint_p: competitive over ligands within the receiver). Checked on
                 shuffled labels: each ligand's number of receivers with p < 0.05 against Binomial(K, 0.05), and
                 how stable a ligand's ranks are between two shuffles (a ligand favoured whatever the labels).
Step 3  filter:  Null B's BH calls whose ligand leaves a footprint in the receiver (p_footprint < 0.05): a subset of
                 Null B's calls, so it adds no call Null B did not make. Also recorded: BH on max(p_null_b,
                 p_footprint) (the intersection-union form; it cannot reject, because p_footprint's floor is
                 1/n_ligands and only one ligand per receiver reaches it) and, as a side result, BH on the Cauchy
                 combination (Liu & Xie 2020) of the two. Judged only on spatial truths (never on NicheNet): TNBC Visium, TNBC Xenium (proxy:
                 ER+/HER2+ tumour, broad classes), Allen brain, hypothalamus MERFISH: odds ratio of colocalisation,
                 Null A's calls kept vs dropped (2,000 bootstrap resamples of (s, t)), as the plan's primary endpoint.

Rows: NATMI from the B = 5,000 addendum (data/liana_nulls_b5000/; at B = 1,000 its Null B makes no BH call on
most data sets), CellChat and CellPhoneDB from benchmark 1's core (B = 1,000; rerun here, rerun_tools.py).
A row whose ligand NicheNet does not have gets p_footprint = 1 (no footprint evidence).
Results: benchmarks/results/liana_nulls/addendum_nichenet.json.
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd
from scipy import stats as sstats

from ccc_nulls import datasets
from ccc_nulls.stats import cauchy_combine

from . import downstream, nichenet
from .common import AMENDMENT_DATASETS, CORE_DATASETS, OUT, log, rows_dir, save
from .core import ALPHA, bh_over

TOOLS = ("natmi", "cellchat", "cellphonedb")
OTHER_NULLS = ("null_b", "b2", "loo_t")      # each judged against Null A's calls (kept vs dropped)
NAME = "addendum_nichenet"
NOTE = "addition, not pre-registered: NicheNet v2 (Zenodo 7074291) as truth (step 1), footprint p (step 2), filter (step 3)"


def rows_path(tool: str, name: str):
    if tool == "natmi":
        return datasets.data_dir() / "liana_nulls_b5000" / f"{name}_natmi.pkl"
    return rows_dir() / f"{name}_{tool}.pkl"


def spatial_truths(name: str) -> dict:
    from .stages import spatial_truth
    out = {}
    if name == "tnbc_wu2021":
        out["visium"] = spatial_truth(name)
        f = OUT / "truth_tnbc_wu2021_xenium.csv"
        if f.exists():
            out["xenium_proxy"] = pd.read_csv(f)[["source", "target", "colocalised"]]
    elif name == "allen_brain_atlas":
        out["visium"] = spatial_truth(name)
    elif name == "moffitt2018_hypothalamus":
        out["merfish"] = pd.read_csv(OUT / f"truth_{name}_spatial.csv")
    return out


def profiles_for(name: str, adata=None, groupby=None, shuffle_seed=None):
    from .data import SPECIES, load, to_human
    if adata is None:
        adata, groupby, _ = load(name, extra_genes=set(nichenet.matrix().index))
        if SPECIES[name] != "human":
            adata, _ = to_human(adata, SPECIES[name])
    if shuffle_seed is not None:
        adata = adata.copy()
        adata.obs[groupby] = pd.Categorical(np.random.default_rng(shuffle_seed).permutation(adata.obs[groupby].astype(str).to_numpy()))
    return downstream.contrast_profiles(adata, groupby), adata, groupby


def calibration(P_shuffles: list) -> dict:
    """Step 2's check on shuffled labels (see module docstring)."""
    out = []
    for P in P_shuffles:
        K = P.shape[1]
        hits = (P.to_numpy() < 0.05).sum(1)
        lim = sstats.binom.ppf(0.99, K, 0.05)
        out.append(dict(n_ligands=int(P.shape[0]), n_receivers=int(K), share_p_lt_005=float((P.to_numpy() < 0.05).mean()),
                        share_ligands_above_binom99=float((hits > lim).mean()), binom99_limit=float(lim)))
    stab = None
    if len(P_shuffles) >= 2:
        a, b = P_shuffles[0].mean(1), P_shuffles[1].mean(1)
        common = a.index.intersection(b.index)
        stab = float(sstats.spearmanr(a[common], b[common]).statistic)
    return dict(shuffles=out, ligand_mean_p_spearman_between_shuffles=stab,
                reading="share_ligands_above_binom99 near 0.01 and a Spearman near 0: no ligand is favoured "
                        "whatever the labels")


def row_footprint_p(rows: pd.DataFrame, P: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    lig = set(P.index)
    L = rows.ligand_complex.map(lambda cx: nichenet.ligand_of(cx, lig)).to_numpy(dtype=object)
    stacked = P.stack()
    keys = list(zip(L, rows.target.astype(str)))
    p = np.array([stacked.get(k, np.nan) if k[0] is not None else np.nan for k in keys], float)
    covered = np.isfinite(p)
    return np.where(covered, p, 1.0), covered


def kept_vs_dropped(rows, truth, call) -> dict:
    from .stages import kept_vs_dropped as kvd
    k = kvd(rows.assign(call_null_b=np.asarray(call, bool)), truth)
    if not (k["kept_calls"] and k["dropped_calls"]):
        k |= dict(odds_ratio=None, ci95=None)
    return k


def evaluate(tool: str, name: str, prof: pd.DataFrame, P: pd.DataFrame) -> dict:
    from .data import SPECIES
    f = rows_path(tool, name)
    if not f.exists():
        return {"missing_rows": str(f)}
    rows = pd.read_pickle(f)
    for c in ("source", "target"):
        rows[c] = rows[c].astype(str)
    if "call_loo_t" not in rows:                              # CellChat, CellPhoneDB: LOO-t from the scores (Amendment 4)
        from .amendment4 import loo_t
        rows["p_loo_t"], rows["call_loo_t"] = loo_t(rows)
    rh = downstream.human_rows(rows, SPECIES[name])          # mouse CellChat rows carry mouse symbols
    res = {"rows": int(len(rows)), "family": int(rows.family.sum()),
           "calls": {n: int(rows[f"call_{n}"].sum()) for n in OTHER_NULLS + ("null_a",) if f"call_{n}" in rows}}

    # step 1: NicheNet as a truth; Null B with D2, B2 and LOO-t with D1 (kept vs dropped among Null A's calls)
    r = downstream.truth_rows(rh, prof, "nichenet")
    res["step1_truth"] = None
    if len(r):
        res["step1_truth"] = {"null_b": downstream.endpoints(r)}
        for n in OTHER_NULLS[1:]:
            if f"call_{n}" in r:
                kept = r[f"call_{n}"].to_numpy(bool)
                res["step1_truth"][n] = {"D1": downstream.odds_ratio(r, kept, ~kept)}

    # step 2: footprint p on the tool's rows
    p2, cov = row_footprint_p(rh, P)
    fam = rows.family.to_numpy(bool)
    a = rows.call_null_a.to_numpy(bool) & (rows.source != rows.target).to_numpy()
    res["step2_footprint"] = dict(covered_rows=int(cov.sum()), covered_family=int((cov & fam).sum()),
                                  share_p_lt_005_null_a_calls=float((p2[a & cov] < 0.05).mean()) if (a & cov).any() else None,
                                  share_p_lt_005_family=float((p2[fam & cov] < 0.05).mean()) if (fam & cov).any() else None)

    # step 3: the filter
    pb = rows.p_null_b.fillna(1.0).to_numpy(float)
    p_and = np.maximum(pb, p2)
    p_or = cauchy_combine(np.column_stack([pb, p2]))
    nb = rows.call_null_b.to_numpy(bool)
    calls = {"null_a": rows.call_null_a.to_numpy(bool), "null_b": nb,
             **{n: rows[f"call_{n}"].to_numpy(bool) for n in OTHER_NULLS[1:] if f"call_{n}" in rows},
             # the filter: Null B's BH calls whose ligand leaves a footprint in the receiver (top 5% of ligands there)
             "filter_null_b_then_footprint": nb & (p2 < ALPHA),
             # documented dead end: BH on max(p_B, p_footprint) cannot reject, since p_footprint's floor is
             # 1/n_ligands and only one ligand per receiver reaches it
             "and_max_bh": bh_over(p_and, fam) <= ALPHA,
             "or_cauchy_bh": bh_over(p_or, fam) <= ALPHA}
    res["step3_calls"] = {k: int(v.sum()) for k, v in calls.items()}
    res["step3_spatial"] = {t: {k: kept_vs_dropped(rows, truth, calls[k]) for k in calls if k != "null_a"}
                            for t, truth in spatial_truths(name).items()}
    return res


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_nulls.nichenet_eval")
    ap.add_argument("--tool", nargs="+", default=list(TOOLS))
    ap.add_argument("--dataset", nargs="+", default=CORE_DATASETS + AMENDMENT_DATASETS)
    ap.add_argument("--shuffles", type=int, default=2)
    args = ap.parse_args(argv)
    f = OUT / f"{NAME}.json"
    out = json.loads(f.read_text()).get("results", {}) if f.exists() else {}
    for name in args.dataset:
        prof, adata, groupby = profiles_for(name)
        P = nichenet.footprint_p(prof)
        d = out.setdefault(name, {})
        if "step2_calibration" not in d:
            d["step2_calibration"] = calibration([nichenet.footprint_p(profiles_for(name, adata, groupby, s)[0])
                                                  for s in range(args.shuffles)])
        for tool in args.tool:
            log(f"nichenet {name} {tool}")
            d[tool] = evaluate(tool, name, prof, P)
            save(NAME, {"addendum": NOTE, "alpha": ALPHA, "results": out})
            r = d[tool]
            if "step3_calls" in r:
                d1 = {n: (v.get("D1") or {}).get("odds_ratio") for n, v in (r["step1_truth"] or {}).items()}
                log(f"  calls {r['step3_calls']}; step 1 D1 {d1}")
        del adata


if __name__ == "__main__":
    main()
