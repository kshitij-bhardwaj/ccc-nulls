"""Our nulls inside LIANA+'s consensus (``rank_aggregate``): the final pipeline's engine.

LIANA+ 1.8.0 ranks interactions by five methods (CellPhoneDB, Connectome, log2FC, NATMI,
SingleCellSignalR) and aggregates their ranks (RRA). Here CellChat is the sixth. CellPhoneDB,
CellChat, and Connectome are the methods with a permutation null, and their p-value (the method's
"specificity") can be replaced by one of ours. Their score (the "magnitude": ``lr_means``,
``lr_probs``, ``expr_prod``) is always LIANA+'s.

    null_a   LIANA+'s own p: all labels permuted together, p = #{s_b >= s} / B (unchanged)
    null_b   the composite null (:mod:`ccc_nulls.composite`) on LIANA+'s score and LIANA+'s own
             permutations: p_L permutes the sender side only (receiver side as observed), p_R the
             receiver side only, p = max(p_L, p_R), with the repository's (1 + k) / (B + 1)
    b2       :mod:`ccc_nulls.competitive` (decoy genes, CAMERA correction) on the same rows

How: LIANA+ scores a permutation method with ``fun(x, perm_stats)``, where ``x`` holds the rows
(source, target, complexes, the observed ligand and receptor statistics) and ``perm_stats`` the
ligand and receptor statistics under each of its label permutations. :func:`capture` and
:func:`capture_connectome` wrap ``fun``: LIANA+'s output is returned untouched, and Null B is
computed from the same ``perm_stats``. One LIANA+ run then gives every variant (:func:`variants`),
each aggregated by LIANA+'s own ``_aggregate``:

    liana5        LIANA+'s rank_aggregate as shipped (five methods)
    liana6        the five plus CellChat, LIANA+'s own p (= swap_null_a)
    swap_<n>      liana6 with CellPhoneDB's, CellChat's, and Connectome's p replaced by null <n>'s
    extra_<n>     liana6 plus methods CellPhoneDB·<n>, CellChat·<n>, Connectome·<n> (the tool's
                  magnitude, null <n>'s p)

A LIANA+ detail inherited here: for a complex, ``perm_stats`` follows the one subunit that is lowest
under the observed labels, as LIANA+'s own permutation test does.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

KEY = ["source", "target", "ligand_complex", "receptor_complex"]
NULLS = ("null_a", "null_b", "b2")
TOOLS = {
    "CellPhoneDB": ("lr_means", "cellphone_pvals"),
    "CellChat": ("lr_probs", "cellchat_pvals"),
    "Connectome": ("expr_prod", "connectome_pvals")
}
CHUNK = 20_000


def _cpdb(L, R):
    """CellPhoneDB's score as LIANA+ computes it for the observed labels."""
    s = (L + R) / 2
    return np.where((L == 0) | (R == 0), 0.0, s)


def _cellchat(L, R):
    from liana.method.sc._cellchat import cellchat
    p = L * R
    return p / (cellchat._kh + p)


def _connectome(L, R):
    """Connectome's score as LIANA+ computes it for the observed labels (expr_prod)."""
    return L * R


SCORES = {
    "CellPhoneDB": (_cpdb, ("ligand_means", "receptor_means")),
    "CellChat": (_cellchat, ("ligand_trimean", "receptor_trimean")),
    "Connectome": (_connectome, ("ligand_means", "receptor_means"))
}


def side_p(score, L, R, perm_stats, chunk: int = CHUNK) -> dict:
    """Null B from LIANA+'s permutations: p_L with the ligand statistic permuted and the receptor as
    observed, p_R the other way round; p = max. ``perm_stats``: (2, B, rows)."""
    L, R = np.asarray(L, float), np.asarray(R, float)
    s = score(L, R)
    B = perm_stats.shape[1]
    kL, kR = np.zeros(len(s), np.int64), np.zeros(len(s), np.int64)
    for a in range(0, len(s), chunk):
        z = slice(a, a + chunk)
        kL[z] = (score(perm_stats[0, :, z], R[None, z]) >= s[None, z]).sum(0)
        kR[z] = (score(L[None, z], perm_stats[1, :, z]) >= s[None, z]).sum(0)
    pval = lambda k: np.where(s == 0, 1.0, (1 + k) / (B + 1))
    pL, pR = pval(kL), pval(kR)
    return dict(score=s, p_L=pL, p_R=pR, p=np.maximum(pL, pR))


def capture(method, store: dict):
    """A copy of LIANA+'s ``method`` (CellPhoneDB or CellChat) whose scoring also writes its rows,
    LIANA+'s p and Null B into ``store[method_name]``. LIANA+'s output is unchanged."""
    from liana.method.sc._Method import Method, MethodMeta

    name = method.method_name
    score, (lcol, rcol) = SCORES[name]

    def fun(x, perm_stats):
        mag, p = method.fun(x, perm_stats)
        nb = side_p(score, x[lcol].to_numpy(), x[rcol].to_numpy(), perm_stats)
        store[name] = pd.DataFrame({**{k: x[k].to_numpy() for k in KEY + ["ligand", "receptor"]},
                                    "magnitude": mag, "p_null_a": p, "p_L_null_b": nb["p_L"],
                                    "p_R_null_b": nb["p_R"], "p_null_b": nb["p"],
                                    "score_check": nb["score"]})
        return mag, p

    meta = MethodMeta(method_name=name, complex_cols=method.complex_cols, add_cols=method.add_cols, fun=fun,
                      magnitude=method.magnitude, magnitude_ascending=method.magnitude_ascending,
                      specificity=method.specificity, specificity_ascending=method.specificity_ascending,
                      permute=method.permute, reference=method.reference)
    return Method(_method=meta)


def capture_connectome(store: dict):
    """A permutation-enabled copy of LIANA+'s Connectome whose scoring writes its rows,
    Null A p and Null B into ``store['Connectome']``."""
    from liana.method.sc._Method import Method, MethodMeta

    score, (lcol, rcol) = SCORES["Connectome"]

    def fun(x, perm_stats):
        L = x[lcol].to_numpy(float)
        R = x[rcol].to_numpy(float)
        mag = score(L, R)
        B = perm_stats.shape[1]
        k_A = (perm_stats[0] * perm_stats[1] >= mag[None, :]).sum(0)
        p_a = np.where(mag == 0, 1.0, (1 + k_A) / (B + 1))
        nb = side_p(score, L, R, perm_stats)
        store["Connectome"] = pd.DataFrame({
            **{k: x[k].to_numpy() for k in KEY + ["ligand", "receptor"]},
            "magnitude": mag,
            "p_null_a": p_a,
            "p_L_null_b": nb["p_L"],
            "p_R_null_b": nb["p_R"],
            "p_null_b": nb["p"],
            "score_check": nb["score"]
        })
        return mag, p_a

    meta = MethodMeta(
        method_name="Connectome",
        complex_cols=["ligand_means", "receptor_means"],
        add_cols=["ligand_zscores", "receptor_zscores"],
        fun=fun,
        magnitude="expr_prod",
        magnitude_ascending=False,
        specificity="connectome_pvals",
        specificity_ascending=True,
        permute=True,
        reference="msraredon/Connectome"
    )
    return Method(_method=meta)


def run(adata, groupby: str, *, resource_name: str = "consensus", n_perms: int = 1000, seed: int = 1337,
        expr_prop: float = 0.1, **kw) -> tuple[dict, dict]:
    """One LIANA+ run of the six methods with CellPhoneDB, Connectome, and CellChat captured. Returns
    ({method name: its table, keys + magnitude + specificity}, {tool: captured rows})."""
    import liana as li
    from liana.method.sc._rank_aggregate import AggregateClass, _rank_aggregate_meta

    store = {}
    methods = [
        capture(li.mt.cellphonedb, store),
        capture_connectome(store),
        li.mt.logfc,
        li.mt.natmi,
        li.mt.singlecellsignalr,
        capture(li.mt.cellchat, store)   # CellChat last: it rescales X in place
    ]
    agg = AggregateClass(_rank_aggregate_meta, methods=methods)
    lrs = agg(adata, groupby=groupby, resource_name=resource_name, n_perms=n_perms, seed=seed, expr_prop=expr_prop,
              consensus_opts=False, inplace=False, use_raw=False, **kw)
    return lrs, store


def add_b2(store: dict, adata, groupby: str, excluded, seed: int = 0) -> None:
    """``p_b2`` (and its sides) for every captured row, from ``adata`` (all genes, for the decoys)."""
    from . import competitive

    types = list(adata.obs[groupby].astype("category").cat.categories.astype(str))
    codes = pd.Categorical(adata.obs[groupby].astype(str), categories=types).codes.astype(np.int64)
    gs = competitive.gene_stats(adata.X, adata.var_names.astype(str), codes, len(types), types=types)
    rows = pd.concat([d[KEY] for d in store.values()]).drop_duplicates().reset_index(drop=True)
    res = competitive.null_b2(adata.X, gs, codes, rows, excluded, seed=seed)
    b2 = pd.concat([rows, res[["gate_b2", "p_L_b2", "p_R_b2", "p_b2"]]], axis=1)
    for name, d in store.items():
        store[name] = d.drop(columns=[c for c in b2.columns if c not in KEY and c in d]).merge(b2, on=KEY, how="left")


def _meta(name, magnitude, specificity):
    from liana.method.sc._Method import MethodMeta
    return MethodMeta(method_name=name, complex_cols=[], add_cols=[], fun=None, magnitude=magnitude,
                      magnitude_ascending=False, specificity=specificity, specificity_ascending=True,
                      permute=False, reference="")


def variants(lrs: dict, store: dict, nulls=None, aggregate_method: str = "rra") -> dict:
    """Every consensus variant, each a LIANA+ consensus table (keys, the methods' scores,
    ``magnitude_rank``, ``specificity_rank``)."""
    import liana as li
    from liana.method._pipe_utils._aggregate import _aggregate
    from liana.method.sc._rank_aggregate import AggregateClass, _rank_aggregate_meta

    if nulls is None:
        # Automatically include b2 only if add_b2 was called
        first_store = next(iter(store.values()))
        nulls = [n for n in NULLS if f"p_{n}" in first_store.columns]

    six = {m.method_name: m for m in [li.mt.cellphonedb, capture_connectome({}), li.mt.logfc, li.mt.natmi,
                                     li.mt.singlecellsignalr, li.mt.cellchat]}
    base = {k: v.copy() for k, v in lrs.items()}

    def agg(tables: dict, metas: list):
        return _aggregate({k: t.copy() for k, t in tables.items()}, consensus=AggregateClass(_rank_aggregate_meta, metas),
                          aggregate_method=aggregate_method).reset_index(drop=True)

    def ours(tool, n):
        mag, spec = TOOLS[tool]
        d = store[tool][KEY + ["magnitude", f"p_{n}"]].drop_duplicates(KEY)
        return d.rename(columns={"magnitude": mag, f"p_{n}": spec})

    out = {"liana5": agg({k: base[k] for k in list(six)[:5]}, list(six.values())[:5]),
           "liana6": agg(base, list(six.values()))}
    for n in nulls:
        if n != "null_a":
            swapped = dict(base)
            for tool in TOOLS:
                swapped[tool] = ours(tool, n)
            out[f"swap_{n}"] = agg(swapped, list(six.values()))
        extra, metas = dict(base), list(six.values())
        for tool, (mag, spec) in TOOLS.items():
            name = f"{tool}·{n}"
            extra[name] = ours(tool, n).rename(columns={mag: f"{mag}__{n}", spec: f"{spec}__{n}"})
            metas.append(_meta(name, f"{mag}__{n}", f"{spec}__{n}"))
        out[f"extra_{n}"] = agg(extra, metas)
    return out
