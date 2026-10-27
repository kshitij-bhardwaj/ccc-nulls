"""Null B against the remaining points of the design checklist (README, "Null B against the checklist").

    python -m tools.cellchat.null_b_checklist --part fake         # TNBC: E2 fake pairs, E4 abundance, A-only by expression
    python -m tools.cellchat.null_b_checklist --part pcal         # Mammary: p-value calibration on known-null rows

Settings are the frozen ones (delta = 0, BH over ``nulls.testable``, B = 1,000); nothing here is tuned.

fake (TNBC, as E2 / E4 of the repository's TNBC benchmark, benchmarks/tnbc_openproblems.py, so Null B
sits in the same table as N1-N5):
- Fake L-R pairs: for every CellChat candidate pair, each ligand and receptor subunit is replaced by
  a background gene (not in CellChatDB, expressed in some cell) from the same mean-expression x
  detection-rate stratum (ccc_nulls.nulls.gene_strata); no cofactors. Complex structure is kept.
  Fake pairs are scored in every (sender, receiver), like the real candidates. They interact in no
  database, so every fake call is false for the claim "s signals to t".
- Unlike the real pairs, fake pairs skip CellChat's over-expression filter (it is run on CellChatDB
  genes only), so they are on average less cell-type specific: the rates compare the two nulls on
  the same rows, not CellChat against the earlier tools.
- One family: BH over the testable real + fake rows. Rates: calls / rows, for real and fake rows.
- E4: calls by decile of the ligand's (and the receptor's) mean expression (min over subunits), real rows.
- A-only (q_A <= 0.05 < q_B) by the same deciles: does Null A's surplus sit at high expression?

pcal (Tabula Sapiens Mammary, Level 1 of the checklist: type-I error and p-value calibration):
- One receptor-shuffle data set (replicate 1 of null_b_calibration.py) and one fully shuffled one.
- Known-null testable rows: P(p <= alpha) for Null A and Null B (a valid p gives <= alpha), and
  quantiles of p for a QQ plot. Calls under Benjamini-Yekutieli (valid under any dependence).

Writes tools/cellchat/results/null_b_checklist_<part>.json.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from ccc_nulls import datasets
from ccc_nulls.nulls import gene_strata
from ccc_nulls.stats import bh
from tools.cellchat import adapter, nulls
from tools.cellchat import cellchat_py as cc
from tools.cellchat.calibrate_null_a import load
from tools.cellchat.null_b_calibration import SEEDS, log, shuffle_genes, side_genes

HERE = Path(__file__).resolve().parent
ALPHA = 0.05
QS = np.linspace(0, 1, 101)


def by(p):
    """Benjamini-Yekutieli q-values: BH times sum_{i<=m} 1/i."""
    return np.minimum(1.0, bh(p) * np.sum(1.0 / np.arange(1, len(p) + 1)))


def adjust(p, rows, fn=bh):
    q = np.ones(len(p))
    if rows.any():
        q[rows] = fn(p[rows])
    return q


def fake_pairs(adata, df, rng):
    """Common-schema rows of fake L-R pairs, one per real candidate pair, in every (sender, receiver)."""
    genes = pd.Index(adata.var_names.astype(str))
    db_genes = set(nulls.signaling_genes())               # every CellChatDB ligand, receptor and cofactor gene
    strata = gene_strata(adata.X)
    expressed = np.asarray((adata.X > 0).sum(0)).ravel() > 0
    pool = expressed & ~genes.isin(list(db_genes))
    pos = {g: i for i, g in enumerate(genes)}
    pairs = df.drop_duplicates("interaction_name")
    st = df[["source", "target"]].drop_duplicates()
    rows = []
    for _, r in pairs.iterrows():
        sub = {}
        for col in ("ligand_genes", "receptor_genes"):
            fake = []
            for g in nulls._genes(r[col]):
                cand = np.where(pool & (strata == strata[pos[g]]))[0]
                fake.append(genes[rng.choice(cand)])
            sub[col] = nulls.SEP.join(fake)
        name = f"FAKE:{r.interaction_name}"
        rows.append(dict(interaction_name=name, ligand=sub["ligand_genes"], receptor=sub["receptor_genes"],
                         ligand_complex=sub["ligand_genes"], receptor_complex=sub["receptor_genes"],
                         ligand_genes=sub["ligand_genes"], receptor_genes=sub["receptor_genes"],
                         agonist_genes="", antagonist_genes="", co_A_genes="", co_I_genes=""))
    fp = pd.DataFrame(rows)
    out = fp.merge(st, how="cross")
    return out[[c for c in df.columns if c in out.columns]]


def decile_of(adata, df, col):
    """Decile (0-9, over real rows' pairs) of the min over subunits of each gene's mean expression."""
    genes = pd.Index(adata.var_names.astype(str))
    mu = np.asarray(adata.X.mean(0)).ravel()
    m = df[col].map(lambda x: min(mu[genes.get_loc(g)] for g in nulls._genes(x))).to_numpy()
    return m


def part_fake(args):
    name = "tnbc_wu2021"
    adata, _ = datasets.load_openproblems(name)
    groupby = datasets.DATASETS[name]["groupby"]
    adata.X = sparse.csr_matrix(adata.X)
    log(f"{name}: {adata.n_obs:,} cells")
    df = adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(adata, groupby), nboot=2, seed=1))
    fk = fake_pairs(adata, df, np.random.default_rng(11))
    both = pd.concat([df.assign(is_fake=False), fk.assign(is_fake=True)], ignore_index=True)
    log(f"real rows {len(df):,} ({df.interaction_name.nunique()} pairs), fake rows {len(fk):,}")
    a = nulls.null_a(adata, both, groupby, B=args.B, seed=0, rng="numpy")
    b = nulls.null_b(adata, both, groupby, B=args.B, seed=0, rng="numpy")
    fake = both.is_fake.to_numpy()
    fam = nulls.testable(a, b)
    q = {"null_a": adjust(a.p_null_a.to_numpy(), fam), "null_b": adjust(b.p_null_b.to_numpy(), fam)}
    out = dict(dataset=name, B=args.B, n_real=int((~fake).sum()), n_fake=int(fake.sum()),
               n_testable_real=int((fam & ~fake).sum()), n_testable_fake=int((fam & fake).sum()), e2={})
    for m, v in q.items():
        c = v <= ALPHA
        out["e2"][m] = dict(real_calls=int((c & ~fake).sum()), fake_calls=int((c & fake).sum()),
                            real_rate=float(c[~fake].mean()), fake_rate=float(c[fake].mean()),
                            real_rate_testable=float(c[fam & ~fake].mean()), fake_rate_testable=float(c[fam & fake].mean()),
                            fake_share_of_calls=float((c & fake).sum() / max(c.sum(), 1)))
    log("E2:", out["e2"])

    real = ~fake
    ca, cb = q["null_a"][real] <= ALPHA, q["null_b"][real] <= ALPHA
    out["e4"] = {}
    for side, col in (("ligand", "ligand_genes"), ("receptor", "receptor_genes")):
        mu = decile_of(adata, both[real], col)
        dec = pd.qcut(pd.Series(mu).rank(method="first"), 10, labels=False).to_numpy()
        out["e4"][side] = [dict(decile=int(d), rows=int((dec == d).sum()),
                                null_a_rate=float(ca[dec == d].mean()), null_b_rate=float(cb[dec == d].mean()),
                                a_only=int((ca & ~cb & (dec == d)).sum()))
                           for d in range(10)]
    log("E4 ligand:", [(r["decile"], round(r["null_a_rate"], 4), round(r["null_b_rate"], 4)) for r in out["e4"]["ligand"]])
    return out


def part_pcal(args):
    adata = load("mammary")
    lig_only, rec_only = side_genes(cc.subset_db(cc.load_db("human")))
    out = dict(dataset="Tabula Sapiens - Mammary", B=args.B, quantiles=QS.tolist())
    for scen in ("receptor", "labels"):
        rng = np.random.default_rng(SEEDS[scen])            # replicate 1 of null_b_calibration.py
        if scen == "labels":
            ad = adata.copy()
            ad.obs["label"] = pd.Categorical(rng.permutation(adata.obs["label"].astype(str).to_numpy()))
        else:
            ad = shuffle_genes(adata, rec_only, rng)
        df = adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(ad, "label"), nboot=2, seed=1))
        a = nulls.null_a(ad, df, "label", B=args.B, seed=1, rng="numpy")
        b = nulls.null_b(ad, df, "label", B=args.B, seed=1, rng="numpy")
        fam = nulls.testable(a, b)
        if scen == "labels":
            known = np.ones(len(df), bool)
        else:
            known = df.receptor_genes.map(lambda x: all(g in rec_only for g in nulls._genes(x))).to_numpy()
        rows = fam & known
        res = dict(n_rows=int(len(df)), n_known_null_testable=int(rows.sum()))
        for m, p in (("null_a", a.p_null_a.to_numpy()), ("null_b", b.p_null_b.to_numpy())):
            pk = p[rows]
            res[m] = dict(type1={str(al): float((pk <= al).mean()) for al in (0.01, 0.05, 0.10)},
                          quantiles=np.quantile(pk, QS).tolist(),
                          bh_false_calls=int(((adjust(p, fam) <= ALPHA) & known).sum()),
                          by_false_calls=int(((adjust(p, fam, by) <= ALPHA) & known).sum()),
                          bh_calls=int((adjust(p, fam) <= ALPHA).sum()), by_calls=int((adjust(p, fam, by) <= ALPHA).sum()))
        out[scen] = res
        log(scen, {m: (res[m]["type1"], res[m]["bh_false_calls"], res[m]["by_false_calls"]) for m in ("null_a", "null_b")})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", required=True, choices=["fake", "pcal"])
    ap.add_argument("--B", type=int, default=1000)
    args = ap.parse_args()
    t0 = time.time()
    out = part_fake(args) if args.part == "fake" else part_pcal(args)
    out["seconds"] = round(time.time() - t0, 1)
    f = HERE / "results" / f"null_b_checklist_{args.part}.json"
    f.write_text(json.dumps(out, indent=2))
    log("wrote", f)


if __name__ == "__main__":
    main()
