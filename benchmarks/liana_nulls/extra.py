"""Amendment 3 (AMENDMENT_3.md): the extra Null B columns, for ``core`` (new runs) and ``extend``
(cached core rows).

    Null B2      ccc_nulls.competitive: decoy genes with the CAMERA correction, both tools.
                 Columns p_L_b2, p_R_b2, p_b2 and the ablation variants p_b2_uncorrected,
                 p_b2_unmatched, p_b2_joint; q_ and call_ over the tool's family.
    Null B (#9)  CellPhoneDB only: benchmarks/null_b_cellphonedb.run unchanged (its side score, its
                 permutations, its three procedures). Columns p_null_b_pr9, q_null_b_pr9 (BH over
                 the gate), q_null_b_pr9_testable, q_null_b_pr9_by and their calls.

B2 reads the same genes the tool read: a mouse data set's CellPhoneDB rows name human orthologs, so
B2 runs on the renamed data there.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccc_nulls import competitive

from .core import ALPHA, KEY, bh_over, keyed

B2_COLUMNS = ["p_b2", "p_b2_uncorrected", "p_b2_unmatched", "p_b2_joint"]
PR9_PROCEDURES = {"null_b_pr9": "null_b", "null_b_pr9_testable": "null_b|testable", "null_b_pr9_by": "null_b|by"}
N_DECOY_UNIVERSE = 5000


def excluded_genes() -> set:
    """Genes of every database used here, both species (decoys must be outside all of them)."""
    from .data import database_genes
    return database_genes("human") | database_genes("mouse")


def b2_columns(adata, groupby: str, rows: pd.DataFrame, excluded: set | None = None, seed: int = 0) -> pd.DataFrame:
    """B2 and its variants for ``rows`` on ``adata`` (X = ln(1 + CP10k)); BH over ``rows.family``."""
    excluded = excluded_genes() if excluded is None else excluded
    types = list(adata.obs[groupby].astype("category").cat.categories.astype(str))
    codes = pd.Categorical(adata.obs[groupby].astype(str), categories=types).codes.astype(np.int64)
    gs = competitive.gene_stats(adata.X, adata.var_names.astype(str), codes, len(types), types=types)
    out = competitive.null_b2(adata.X, gs, codes, rows, excluded, seed=seed)
    fam = rows.family.to_numpy(bool)
    for c in ["p_b2"] + B2_COLUMNS[1:]:
        name = c[2:]
        out[f"q_{name}"] = bh_over(out[c].to_numpy(), fam)
        out[f"call_{name}"] = out[f"q_{name}"].to_numpy() <= ALPHA
    return out


def pr9_columns(adata, groupby: str, rows: pd.DataFrame, workdir: Path, B: int = 1000, seed: int = 0) -> tuple[pd.DataFrame, dict]:
    """#9's CellPhoneDB Null B, joined onto ``rows`` by (source, target, ligand, receptor complex)."""
    from benchmarks import null_b_cellphonedb as pr9

    from .common import cpdb_database

    res, q, meta, _ = pr9.run(adata, groupby, str(cpdb_database()), B, seed, Path(workdir))
    res = res.assign(p_null_b_pr9=res.p_null_b.to_numpy(),
                     **{f"q_{k}": q[v] for k, v in PR9_PROCEDURES.items()})
    res = keyed(res)
    cols = ["p_null_b_pr9"] + [f"q_{k}" for k in PR9_PROCEDURES]
    j = keyed(rows[KEY]).merge(res[KEY + cols].drop_duplicates(KEY), on=KEY, how="left")
    out = pd.DataFrame(index=rows.index)
    for c in cols:
        out[c] = j[c].fillna(1.0).to_numpy()
    for k in PR9_PROCEDURES:
        out[f"call_{k}"] = out[f"q_{k}"].to_numpy() <= ALPHA
    meta = {k: v for k, v in meta.items()} | dict(rows_matched=int(j.p_null_b_pr9.notna().sum()), rows=int(len(rows)), B=B, seed=seed)
    return out, meta


def agreement(a, b) -> dict:
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    u = (a | b).sum()
    return dict(calls_a=int(a.sum()), calls_b=int(b.sum()), same_call=float((a == b).mean()),
                jaccard=float((a & b).sum() / u) if u else float("nan"))


def add_extra(adata, groupby: str, tool, df: pd.DataFrame, workdir: Path, species: str, log=print) -> tuple[pd.DataFrame, dict]:
    """The Amendment 3 columns for one tool's rows (``df`` already has ``family``, ``call_null_b``)."""
    from scipy import stats as sstats

    from .data import to_human

    meta = {}
    a = adata
    if tool.name == "cellphonedb" and species != "human":
        a, _ = to_human(adata, species)
    b2 = b2_columns(a, groupby, df)
    for c in b2.columns:
        df[c] = b2[c].to_numpy()
    meta["b2"] = dict(gated_rows=int(df.gate_b2.sum()), gated_in_family=int((df.gate_b2 & df.family).sum()), family=int(df.family.sum()),
                      calls={v: int(df[f"call_{v[2:]}"].sum()) for v in B2_COLUMNS},
                      vs_null_b=agreement(df.call_b2, df.call_null_b), vs_null_a=agreement(df.call_b2, df.call_null_a))
    log(f"    B2 calls {meta['b2']['calls']}")
    if tool.name == "cellphonedb":
        pr, m = pr9_columns(a, groupby, df, Path(workdir) / "pr9")
        for c in pr.columns:
            df[c] = pr[c].to_numpy()
        fam = df.family.to_numpy(bool)
        rho = sstats.spearmanr(df.p_null_b_pr9[fam], df.p_null_b[fam]).statistic if fam.any() else float("nan")
        meta["pr9"] = m | dict(calls={k: int(df[f"call_{k}"].sum()) for k in PR9_PROCEDURES},
                               vs_null_b=agreement(df.call_null_b_pr9, df.call_null_b),
                               p_identical_in_family=float((df.p_null_b_pr9[fam] == df.p_null_b[fam]).mean()),
                               spearman_p_in_family=float(rho))
        log(f"    #9 calls {meta['pr9']['calls']}; vs ours {meta['pr9']['vs_null_b']}")
    if a is not adata:
        del a
    return df, meta
