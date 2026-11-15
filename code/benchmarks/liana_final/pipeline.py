"""One data set through the final pipeline: LIANA+'s consensus with our nulls (PLAN.md).

``consensus(adata, groupby, species)`` runs LIANA+ once (six methods, CellPhoneDB and CellChat
captured; :mod:`ccc_nulls.liana_consensus`), adds B2 and returns every consensus variant. ``wide``
puts the variants side by side, ``fidelity`` is gate F0.
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from ccc_nulls import datasets
from ccc_nulls import liana_consensus as lc

from benchmarks.liana_nulls.common import log

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results" / "liana_final"
ALPHA = 0.05
N_PERMS, SEED = 1000, 1337
VARIANTS = ["liana5", "liana6", "swap_null_b", "swap_b2", "extra_null_a", "extra_null_b", "extra_b2"]
JUDGED = ["swap_null_b", "extra_null_b"]
BASELINE = "liana6"


def rows_dir() -> Path:
    d = Path(datasets.data_dir()) / "liana_final"
    d.mkdir(parents=True, exist_ok=True)
    return d


def plan_hash() -> str:
    return hashlib.sha256((HERE / "PLAN.md").read_bytes()).hexdigest()


def amendment_hash(n: int = 1) -> str | None:
    f = HERE / f"AMENDMENT_{n}.md"
    return hashlib.sha256(f.read_bytes()).hexdigest() if f.exists() else None


def save(name: str, obj: dict):
    OUT.mkdir(parents=True, exist_ok=True)
    head = {"plan_sha256": plan_hash()} | ({"amendment_1_sha256": h} if (h := amendment_hash(1)) else {})
    obj = {**head, "written": time.strftime("%Y-%m-%d %H:%M"), **obj}
    (OUT / f"{name}.json").write_text(json.dumps(obj, indent=1, default=str))
    log("wrote", OUT / f"{name}.json")


def resource(species: str) -> str:
    return "consensus" if species == "human" else f"{species}consensus"


def liana_view(adata, extra_genes=()):
    """The cells with the genes LIANA+ saw in benchmarks/liana_nulls (Tabula Sapiens' decoy genes are
    for B2 only: CellChat's scaling reads the maximum of X)."""
    if not len(extra_genes):
        return adata
    return adata[:, ~adata.var_names.isin(list(extra_genes))].copy()


def consensus(adata, groupby: str, species: str, extra_genes=(), n_perms: int = N_PERMS, seed: int = SEED):
    """({variant: LIANA+ consensus table}, {tool: captured rows with every null's p}, meta)."""
    from benchmarks.liana_nulls.extra import excluded_genes

    t0 = time.time()
    a = liana_view(adata, extra_genes)
    lrs, store = lc.run(a.copy(), groupby, resource_name=resource(species), n_perms=n_perms, seed=seed)
    t1 = time.time()
    lc.add_b2(store, adata, groupby, excluded_genes())
    V = lc.variants(lrs, store)
    meta = dict(resource=resource(species), n_perms=n_perms, seed=seed, n_cells=int(adata.n_obs),
                n_types=int(adata.obs[groupby].nunique()), seconds_liana=round(t1 - t0, 1),
                seconds_total=round(time.time() - t0, 1), rows={k: int(len(v)) for k, v in V.items()})
    return V, store, meta


def wide(V: dict) -> pd.DataFrame:
    """KEY + magnitude_rank__v, specificity_rank__v, call__v for every variant."""
    out = None
    for v, t in V.items():
        t = t[lc.KEY + ["magnitude_rank", "specificity_rank"]].drop_duplicates(lc.KEY)
        t = t.rename(columns={"magnitude_rank": f"magnitude_rank__{v}", "specificity_rank": f"specificity_rank__{v}"})
        out = t if out is None else out.merge(t, on=lc.KEY, how="outer")
    for v in V:
        out[f"call__{v}"] = out[f"specificity_rank__{v}"].to_numpy() <= ALPHA
    for c in ("source", "target"):
        out[c] = out[c].astype(str)
    return out


def _same(ours: pd.DataFrame, ref: pd.DataFrame, cols: list) -> dict:
    m = ours.drop_duplicates(lc.KEY).merge(ref.drop_duplicates(lc.KEY), on=lc.KEY, suffixes=("_ours", "_liana"))
    return dict(rows_ours=int(len(ours)), rows_liana=int(len(ref)), rows_matched=int(len(m)),
                identical={c: bool(np.array_equal(m[f"{c}_ours"].to_numpy(), m[f"{c}_liana"].to_numpy())) for c in cols})


def fidelity(adata, groupby: str, species: str, V: dict, store: dict, extra_genes=(), standalone: bool = False,
             n_perms: int = N_PERMS, seed: int = SEED) -> dict:
    """F0: liana5 against LIANA+'s rank_aggregate; with ``standalone``, the captured tools against
    LIANA+'s own CellPhoneDB and CellChat."""
    import liana as li

    a = liana_view(adata, extra_genes)
    kw = dict(resource_name=resource(species), n_perms=n_perms, seed=seed, use_raw=False, inplace=False)
    ref = li.mt.rank_aggregate(a.copy(), groupby, **kw)
    out = {"liana5_vs_rank_aggregate": _same(V["liana5"], ref, ["magnitude_rank", "specificity_rank", "lr_means",
                                                                 "cellphone_pvals", "expr_prod", "lrscore"])}
    if standalone:
        for tool, m in (("CellPhoneDB", li.mt.cellphonedb), ("CellChat", li.mt.cellchat)):
            mag, spec = lc.TOOLS[tool]
            s = m(a.copy(), groupby, **kw)
            out[f"{tool}_vs_standalone"] = _same(store[tool].rename(columns={"magnitude": mag, "p_null_a": spec}), s, [mag, spec])
    out["passed"] = all(all(v["identical"].values()) and v["rows_ours"] == v["rows_liana"] == v["rows_matched"]
                        for k, v in out.items() if isinstance(v, dict))
    return out


def summary(tbl: pd.DataFrame, store: dict) -> dict:
    """Calls per variant, overlap with the baseline, and how our p compares with LIANA+'s."""
    base = tbl[f"call__{BASELINE}"].to_numpy(bool)
    calls = {}
    for v in VARIANTS:
        c = tbl[f"call__{v}"].to_numpy(bool)
        u = (c | base).sum()
        calls[v] = dict(calls=int(c.sum()), jaccard_with_liana6=float((c & base).sum() / u) if u else float("nan"),
                        kept_of_liana6=int((c & base).sum()), outside_liana6=int((c & ~base).sum()))
    tools = {}
    for tool, d in store.items():
        tools[tool] = dict(rows=int(len(d)), p_null_a_lt_05=int((d.p_null_a < ALPHA).sum()),
                           p_null_b_lt_05=int((d.p_null_b < ALPHA).sum()), p_b2_lt_05=int((d.p_b2 < ALPHA).sum()),
                           b2_gate=float(d.gate_b2.mean()), null_b_ge_null_a=float((d.p_null_b >= d.p_null_a).mean()),
                           score_max_abs_diff=float(np.abs(d.score_check - d.magnitude).max()))
    return dict(rows=int(len(tbl)), calls=calls, tools=tools)
