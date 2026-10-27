"""Null B calibration on Tabula Sapiens: known-null data for the claim Null B tests. TNBC is not touched.

    python -m tools.cellchat.null_b_calibration --tissue mammary --scenario receptor --reps 20
    python -m tools.cellchat.null_b_calibration --summarise

Design (each choice and why):
- The claim is "s signals to t": the sender side is specific to s AND the receiver side to t. Known-null
  data must make that claim false while keeping everything else, so three scenarios, each replicated
  (20 independent null data sets per tissue, as in calibrate_null_a.py):
    labels    every cell's label shuffled: both sides null (the standard check; also Null A's null).
    receptor  the receptor-side genes (receptors, co-receptors) shuffled across cells with one cell
              permutation: their joint distribution (complexes) is kept, their link to cell type is
              removed, the ligand side is untouched. Every L-R pair whose receiver-side genes are all
              shuffled is null for the claim, in every (sender, receiver).
    ligand    the mirror image: ligand-side genes shuffled.
  Genes on both sides (a gene that is a ligand in one pair and a receptor in another, agonists and
  antagonists, which enter both sides) are never shuffled, and pairs that use them are not "known
  null"; they stay in the data and in the multiple-testing family.
- Same pipeline, candidates, observed prob and permutations for both nulls: CellChat's pipeline once
  per data set, then null_a and null_b with the same seed (B = 1,000).
- Decision: Benjamini-Hochberg at alpha in {0.01, 0.05, 0.10} over the testable rows
  (``nulls.testable``: some labelling gives prob > 0), whose membership does not depend on which
  labelling is the observed one. Ablation: BH over the rows CellChat reports (prob > 0), the family
  users correct over, which selects on the observed labels.
- Error rates. On "labels" every call is false: FDR = P(at least one call), Clopper-Pearson 95% CI.
  On "receptor" / "ligand", V = calls among known-null rows (all false), R = all calls. The FDR is
  E[V / R] (Benjamini & Hochberg 1995): its estimate is the mean over replicates of V / R, with a
  t-interval; it is a lower bound, since calls on the other rows may be false too. P(V > 0), the
  chance of any false call, is reported as well but is not what FDR control bounds: with real
  signal in the other rows, BH at alpha is allowed about alpha x R false calls.
- Procedures: Null A (CellChat's test, BH), Null B (max(p_L, p_R), BH), Null B at delta = 0.01 and
  0.05 (interval null, sensitivity), and each side alone (p_L or p_R: ablations showing that both
  halves are needed).
- Real labels: calls of each procedure; calls made by Null A but not Null B, by the side that fails.

Writes tools/cellchat/results/null_b_calibration_<tissue>_<scenario>.json and, with --summarise,
null_b_calibration_summary.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats as sstats

from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat import cellchat_py as cc
from tools.cellchat.calibrate_null_a import TISSUES, load

HERE = Path(__file__).resolve().parent
ALPHAS = (0.01, 0.05, 0.10)
DELTAS = (0.0, 0.01, 0.05)
SEEDS = {"labels": 1000, "receptor": 2000, "ligand": 3000}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def side_genes(db) -> tuple[set, set]:
    """CellChatDB genes that are only ever on the sender side (ligands) or only on the receiver side
    (receptors, co-receptors). Agonists / antagonists enter both sides and are on neither list."""
    inter, cpx, cof, info = db["interaction"], db["complex"], db["cofactor"], db["geneInfo"]

    def cofactors(cols):
        names = dict.fromkeys(c for col in cols for c in inter[col] if c)
        sub = cof.loc[[c for c in names if c in cof.index]].filter(like="cofactor")
        return {g for g in sub.to_numpy().ravel() if g}

    lig = set(cc._extract_gene_subset(list(dict.fromkeys(inter.ligand)), cpx, info))
    rec = set(cc._extract_gene_subset(list(dict.fromkeys(inter.receptor)), cpx, info))
    rec |= cofactors(["co_A_receptor", "co_I_receptor"])
    both = cofactors(["agonist", "antagonist"])
    return lig - rec - both, rec - lig - both


def shuffle_genes(adata, genes: set, rng):
    """Copy of ``adata`` with the columns of ``genes`` moved to other cells by one random permutation."""
    m = adata.var_names.isin(list(genes)).astype(float)
    X = sparse.csr_matrix(adata.X)
    out = adata.copy()
    out.X = sparse.csr_matrix(X @ sparse.diags(1 - m) + X[rng.permutation(adata.n_obs)] @ sparse.diags(m))
    return out


def one_dataset(adata, B, seed, known_null_side=None, shuffled=frozenset()):
    """Pipeline, Null A and Null B on the same permutations; calls of every procedure."""
    df = adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(adata, "label"), nboot=2, seed=1))
    t0 = time.time()
    a = nulls.null_a(adata, df, "label", B=B, seed=seed, rng="numpy")
    b = nulls.null_b(adata, df, "label", B=B, seed=seed, rng="numpy", delta=DELTAS)
    secs = time.time() - t0
    assert np.array_equal(a.score_null_a.to_numpy(), b.score_null_b.to_numpy())
    scored = b.score_null_b.to_numpy() > 0
    family = nulls.testable(a, b)
    p = {"null_a": a.p_null_a.to_numpy(), "null_b": b.p_null_b.to_numpy(),
         **{f"null_b|delta={d:g}": b[f"p_null_b_delta={d:g}"].to_numpy() for d in DELTAS[1:]},
         "sender_side_only": b.p_L_null_b.to_numpy(), "receiver_side_only": b.p_R_null_b.to_numpy()}
    q = {}
    for m, v in [*p.items(), ("null_a|reported_rows", p["null_a"]), ("null_b|reported_rows", p["null_b"])]:
        rows = scored if m.endswith("reported_rows") else family
        q[m] = np.ones(len(df))
        if rows.any():
            q[m][rows] = bh(v[rows])
    if known_null_side is None:
        null_rows = np.zeros(len(df), bool)
    else:
        col = "receptor_genes" if known_null_side == "receptor" else "ligand_genes"
        null_rows = df[col].map(lambda x: all(g in shuffled for g in nulls._genes(x))).to_numpy()
    out = dict(n_lr=int(df.interaction_name.nunique()), n_rows=int(len(df)), n_scored=int(scored.sum()),
               n_testable=int(family.sum()), n_known_null_scored=int((null_rows & scored).sum()), seconds=round(secs, 1), calls={})
    for m in q:
        out["calls"][m] = {str(al): dict(R=int((q[m] <= al).sum()), V=int(((q[m] <= al) & null_rows).sum()))
                           for al in ALPHAS}
    called_a, called_b = q["null_a"] <= 0.05, q["null_b"] <= 0.05
    side = b.limiting_side_null_b.to_numpy()
    only_a = called_a & ~called_b
    out["a_only"] = dict(n=int(only_a.sum()), sender_side_fails=int((only_a & (side == "L")).sum()),
                         receiver_side_fails=int((only_a & (side == "R")).sum()),
                         median_prob=float(np.median(b.score_null_b.to_numpy()[only_a])) if only_a.any() else None)
    out["b_only"] = int((called_b & ~called_a).sum())
    out["both"] = int((called_a & called_b).sum())
    if called_b.any():
        out["median_prob_both"] = float(np.median(b.score_null_b.to_numpy()[called_a & called_b]))
    return out


def clopper_pearson(k, n):
    lo = sstats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = sstats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return [float(lo), float(hi)]


def summarise_reps(reps, scenario):
    """Per procedure and alpha. labels: FDR = P(R > 0). receptor / ligand: P(V > 0), mean V / R."""
    out = {}
    for m in reps[0]["calls"]:
        out[m] = {}
        for al in map(str, ALPHAS):
            R = np.array([r["calls"][m][al]["R"] for r in reps])
            V = np.array([r["calls"][m][al]["V"] for r in reps]) if scenario != "labels" else R
            k = int((V > 0).sum())
            fdp = np.where(R > 0, V / np.maximum(R, 1), 0.0)
            half = sstats.t.ppf(0.975, len(fdp) - 1) * fdp.std(ddof=1) / np.sqrt(len(fdp)) if len(fdp) > 1 else np.nan
            out[m][al] = dict(fdr=float(fdp.mean()), fdr_ci95=[float(fdp.mean() - half), float(fdp.mean() + half)],
                              p_any_false_call=k / len(reps), p_any_false_call_ci95=clopper_pearson(k, len(reps)),
                              mean_false_calls=float(V.mean()), mean_calls=float(R.mean()), n_datasets=len(reps))
    return out


def summarise_all():
    res = {}
    for f in sorted((HERE / "results").glob("null_b_calibration_*_*.json")):
        if f.stem.endswith("summary"):
            continue
        d = json.loads(f.read_text())
        res.setdefault(d["tissue"], {})[d["scenario"] + d.get("tag", "")] = d.get("summary", d.get("real"))
    out = HERE / "results" / "null_b_calibration_summary.json"
    out.write_text(json.dumps(res, indent=2))
    log("wrote", out)
    for t, sc in res.items():
        for s, d in sc.items():
            if s.startswith("real"):
                log(t, s, {m: v["0.05"]["R"] for m, v in d["calls"].items()}, "A only:", d["a_only"])
            else:
                log(t, s, {m: (round(v["0.05"]["fdr"], 4), round(v["0.05"]["mean_calls"], 1)) for m, v in d.items()})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", choices=list(TISSUES))
    ap.add_argument("--scenario", choices=["real", *SEEDS])
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--tag", default="")
    ap.add_argument("--summarise", action="store_true")
    args = ap.parse_args()
    if args.summarise:
        return summarise_all()
    adata = load(args.tissue)
    lig_only, rec_only = side_genes(cc.subset_db(cc.load_db("human")))
    f = HERE / "results" / f"null_b_calibration_{args.tissue}_{args.scenario}{args.tag}.json"
    out = dict(dataset=TISSUES[args.tissue]["title"], tissue=args.tissue, scenario=args.scenario, tag=args.tag,
               B=args.B, alphas=ALPHAS, deltas=DELTAS, n_cells=int(adata.n_obs),
               n_ligand_only_genes=len(lig_only & set(adata.var_names)),
               n_receptor_only_genes=len(rec_only & set(adata.var_names)))
    log(f"{args.tissue} / {args.scenario}: {adata.n_obs:,} cells")
    if args.scenario == "real":
        out["real"] = one_dataset(adata, args.B, seed=0)
        log("calls at 0.05:", {m: v["0.05"]["R"] for m, v in out["real"]["calls"].items()})
    else:
        reps = []
        for r in range(args.reps):
            rng = np.random.default_rng(SEEDS[args.scenario] + r)
            if args.scenario == "labels":
                ad = adata.copy()
                ad.obs["label"] = pd.Categorical(rng.permutation(adata.obs.label.to_numpy()))
                reps.append(one_dataset(ad, args.B, seed=r + 1))
            else:
                genes = rec_only if args.scenario == "receptor" else lig_only
                reps.append(one_dataset(shuffle_genes(adata, genes, rng), args.B, seed=r + 1,
                                        known_null_side=args.scenario, shuffled=frozenset(genes)))
            c = reps[-1]["calls"]
            log(f"  rep {r + 1}/{args.reps}: false calls at 0.05: Null A {c['null_a']['0.05']['V' if args.scenario != 'labels' else 'R']}, "
                f"Null B {c['null_b']['0.05']['V' if args.scenario != 'labels' else 'R']} "
                f"(all calls A {c['null_a']['0.05']['R']}, B {c['null_b']['0.05']['R']})")
            out["reps"] = reps
            out["summary"] = summarise_reps(reps, args.scenario)
            f.write_text(json.dumps(out, indent=2))                      # partial results survive a stop
    f.write_text(json.dumps(out, indent=2))
    log("wrote", f)


if __name__ == "__main__":
    main()
