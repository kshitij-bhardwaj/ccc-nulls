"""LIANA+ as it is, plus Null A (native and ours) and Null B for each tool, written back into liana_res.

``run_dataset(adata, groupby, tools)``:
1. ``li.mt.rank_aggregate`` (LIANA+'s consensus, unchanged) -> ``adata.uns["liana_res"]``.
2. Each tool's LIANA+ re-implementation (the extra column).
3. ``Tool.run``: columns 1-3 (native, our Null A, Null B); BH over the tool's family for both nulls.
4. Amendment 3 (``extra.add_extra``): Null B2 for both tools and #9's Null B for CellPhoneDB.
5. Our columns joined onto ``liana_res`` and the tool's LIANA+ table by (source, target, ligand
   complex, receptor complex); the join coverage is reported, never assumed.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from ccc_nulls import liana_reference
from ccc_nulls.stats import bh

from .tools import complex_key

ALPHA = 0.05
KEY = ["source", "target", "ligand_complex", "receptor_complex"]
OUR_COLUMNS = ["p_native", "call_native", "p_port_native", "p_null_a", "q_null_a", "call_null_a",
               "p_null_b", "q_null_b", "call_null_b", "limiting_partner", "family",
               "p_b2", "q_b2", "call_b2",                                        # Amendment 3
               "p_null_b_pr9", "q_null_b_pr9", "call_null_b_pr9"]


def bh_over(p, family):
    q = np.ones(len(p))
    q[family] = bh(np.asarray(p, float)[family])
    return q


def add_q(df: pd.DataFrame, alpha: float = ALPHA) -> pd.DataFrame:
    fam = df.family.to_numpy(bool)
    for null in ("null_a", "null_b"):
        df[f"q_{null}"] = bh_over(df[f"p_{null}"].to_numpy(), fam)
        df[f"call_{null}"] = df[f"q_{null}"].to_numpy() <= alpha
    return df


def keyed(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in ("source", "target"):
        out[c] = out[c].astype(str)
    out["ligand_complex"], out["receptor_complex"] = complex_key(out.ligand_complex), complex_key(out.receptor_complex)
    return out


def liana_builtin(adata, groupby: str, tool, n_perms: int = 1000, seed: int = 0) -> pd.DataFrame:
    """The tool as LIANA+ re-implements it, on the tool's own resource."""
    ref = liana_reference.run(adata, tool.liana_method, groupby, n_perms=n_perms, seed=seed,
                              resource_name=tool.liana_resource)
    ref = keyed(ref)
    ref = ref.drop_duplicates(KEY)
    ref["p_liana"] = ref.p_native
    ref["call_liana"] = ref.p_liana.to_numpy() < ALPHA
    ref["q_liana"] = bh(ref.p_liana.fillna(1.0).to_numpy())
    return ref[KEY + ["score", "p_liana", "call_liana", "q_liana"]].rename(columns={"score": "score_liana"})


def jaccard(a, b) -> float:
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    u = (a | b).sum()
    return float((a & b).sum() / u) if u else float("nan")


def summary(df: pd.DataFrame, tool, meta: dict) -> dict:
    """The three columns (plus LIANA+'s), as numbers of calls and how they overlap."""
    fam = df.family.to_numpy(bool)
    out = {
        "tool": tool.name, "native_rule": tool.native_rule, "rows": int(len(df)), "family": int(fam.sum()),
        "calls": {
            "1_native_tool_rule": int(df.call_native.sum()) if df.p_native.notna().any() else None,
            "2_our_null_a_tool_rule": int(df.call_port_native.sum()),
            "2_our_null_a_bh": int(df.call_null_a.sum()),
            "3_null_b_bh": int(df.call_null_b.sum()),
        },
        "null_b_calls_outside_null_a": int((df.call_null_b & ~df.call_null_a).sum()),
        "null_a_only_by_limiting_partner": df.loc[df.call_null_a & ~df.call_null_b, "limiting_partner"].value_counts().to_dict(),
        "fidelity": meta.get("fidelity"),
        "meta": {k: v for k, v in meta.items() if k != "fidelity"},
    }
    if "p_liana" in df:
        m = df.p_liana.notna().to_numpy()
        out["liana_builtin"] = {
            "rows_matched": int(m.sum()),
            "calls_p_lt_0.05": int(df.call_liana.fillna(False).sum()),
            "jaccard_with_native_rule": jaccard(df.call_liana.fillna(False) & m, df.call_port_native & m),
            "spearman_score": float(pd.Series(df.score[m]).corr(pd.Series(df.score_liana[m]), method="spearman")) if m.any() else None,
        }
    return out


def write_back(liana_res: pd.DataFrame, tool_name: str, df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Our columns for one tool, prefixed with the tool's name, joined onto LIANA+'s consensus table."""
    lr = keyed(liana_res)
    cols = [c for c in OUR_COLUMNS if c in df]
    tmp = [f"_key_{k}" for k in KEY]
    ours = df[KEY + cols].drop_duplicates(KEY).rename(columns={c: f"{tool_name}_{c}" for c in cols})
    ours.columns = tmp + list(ours.columns[len(KEY):])
    old = [c for c in liana_res.columns if c.startswith(f"{tool_name}_") and c[len(tool_name) + 1:] in OUR_COLUMNS]
    merged = liana_res.drop(columns=old).reset_index(drop=True).assign(**{t: lr[k].to_numpy() for t, k in zip(tmp, KEY)})
    merged = merged.merge(ours, on=tmp, how="left").drop(columns=tmp)
    hit = merged[f"{tool_name}_p_null_a"].notna()
    in_liana = ours.merge(lr[KEY].drop_duplicates().set_axis(tmp, axis=1), on=tmp)
    cov = dict(liana_rows=int(len(liana_res)), with_our_columns=int(hit.sum()),
               our_family_rows_in_liana=int(in_liana[f"{tool_name}_family"].sum()),
               our_family_rows=int(df.family.sum()))
    return merged, cov


def run_dataset(adata, groupby: str, tools: list, workdir: Path, *, liana_aggregate: bool = True,
                liana_builtin_column: bool = True, liana_resource: str = "consensus", extra_nulls: bool = True,
                extra_adata=None, log=print) -> tuple[dict, dict]:
    """Returns ({tool: rows}, summary dict); ``adata.uns["liana_res"]`` gains our columns.
    ``extra_adata``: the same cells with more genes (Tabula Sapiens' decoy universe), for Amendment 3
    only, so LIANA+ and the tools see exactly the genes they saw before."""
    import liana as li

    out, rows = {"n_cells": int(adata.n_obs), "n_types": int(adata.obs[groupby].nunique())}, {}
    if liana_aggregate:
        t0 = time.time()
        li.mt.rank_aggregate(adata, groupby=groupby, resource_name=liana_resource, use_raw=False, n_perms=1000,
                             seed=0, verbose=False)
        out["liana_rank_aggregate"] = dict(rows=int(len(adata.uns["liana_res"])), resource=liana_resource,
                                           seconds=round(time.time() - t0, 1))
        log(f"  LIANA+ rank_aggregate: {len(adata.uns['liana_res']):,} rows")
    for tool in tools:
        log(f"  {tool.name}")
        df, meta = tool.run(adata, groupby, Path(workdir) / tool.name, log=log)
        df = add_q(keyed(df))
        if liana_builtin_column:
            t0 = time.time()
            ref = liana_builtin(adata, groupby, tool)
            df = df.merge(ref, on=KEY, how="left")
            meta["seconds_liana_builtin"] = round(time.time() - t0, 1)
        if extra_nulls:
            from .extra import add_extra
            df, meta["amendment_3"] = add_extra(adata if extra_adata is None else extra_adata, groupby, tool, df, Path(workdir) / tool.name, tool.species, log=log)
        out[tool.name] = summary(df, tool, meta)
        if liana_aggregate:
            adata.uns["liana_res"], out[tool.name]["write_back"] = write_back(adata.uns["liana_res"], tool.name, df)
        rows[tool.name] = df
        log(f"    calls: {out[tool.name]['calls']}; fidelity: {meta.get('fidelity')}")
    return rows, out
