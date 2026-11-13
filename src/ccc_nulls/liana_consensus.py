"""Our nulls inside LIANA+'s consensus (``rank_aggregate``): the final pipeline's engine.

LIANA+ 1.8.0 ranks interactions by five methods (CellPhoneDB, Connectome, log2FC, NATMI,
SingleCellSignalR) and aggregates their ranks (RRA). Here CellChat is the sixth. CellPhoneDB and
CellChat are the two with a permutation null, and their p-value (the method's "specificity") can be
replaced by one of ours. Their score (the "magnitude": ``lr_means``, ``lr_probs``) is always LIANA+'s. The nulls::

    null_a   LIANA+'s own p: all labels permuted together, p = #{s_b >= s} / B (unchanged)
    null_b   the composite null (:mod:`ccc_nulls.composite`) on LIANA+'s score and LIANA+'s own
             permutations: p_L permutes the sender side only (receiver side as observed), p_R the
             receiver side only, p = max(p_L, p_R), with the repository's (1 + k) / (B + 1)
    b2       :mod:`ccc_nulls.competitive` (decoy genes, CAMERA correction) on the same rows

How: LIANA+ scores a permutation method with ``fun(x, perm_stats)``, where ``x`` holds the rows
(source, target, complexes, the observed ligand and receptor statistics) and ``perm_stats`` the
ligand and receptor statistics under each of its label permutations. :func:`capture` wraps ``fun``:
LIANA+'s output is returned untouched, and Null B is computed from the same ``perm_stats``. One
LIANA+ run then gives every variant (:func:`variants`), each aggregated by LIANA+'s own ``_aggregate``::

    liana5        LIANA+'s rank_aggregate as shipped (five methods)
    liana6        the five plus CellChat, LIANA+'s own p (= swap_null_a)
    swap_<n>      liana6 with CellPhoneDB's and CellChat's p replaced by null <n>'s
    extra_<n>     liana6 plus two methods, CellPhoneDB·<n> and CellChat·<n> (the tool's magnitude,
                  null <n>'s p), as the 7th and 8th
    swap3_<n>     with ``run(..., connectome=True)``: as swap_<n>, and Connectome's specificity (LIANA+'s
                  ``scaled_weight``, which is not a p-value) replaced by null <n>'s p on Connectome's
                  score; swap3_null_a gives Connectome a p from all labels permuted together

Connectome (idea from PR #13): LIANA+ ranks it by ``expr_prod`` and ``scaled_weight`` and runs no
permutations for it. :func:`capture_connectome` lets LIANA+ permute for it as for CellPhoneDB (same seed,
same means), returns LIANA+'s own two scores unchanged, and stores the nulls' p on ``expr_prod``. So
liana5, liana6 and every ``swap_<n>`` and ``extra_<n>`` variant are the same with or without it.

A LIANA+ detail inherited here: for a complex, ``perm_stats`` follows the one subunit that is lowest
under the observed labels, as LIANA+'s own permutation test does.

``strata`` (e.g. donor per cell) restricts every permutation to within a stratum
(:func:`stratified_permutations`), for Null A and Null B alike. :func:`rank_aggregate` is the
user-facing call: LIANA+'s ``rank_aggregate`` with our null in place of the tools'.
"""
from __future__ import annotations

from contextlib import contextmanager

import numpy as np
import pandas as pd

KEY = ["source", "target", "ligand_complex", "receptor_complex"]
NULLS = ("null_a", "null_b", "b2")
TOOLS = {"CellPhoneDB": ("lr_means", "cellphone_pvals"), "CellChat": ("lr_probs", "cellchat_pvals")}
CONNECTOME = ("expr_prod", "connectome_p")      # Connectome's magnitude, and the column our p goes to
CHUNK = 20_000
ROWS_AT_ONCE = 100_000      # rows scored at once: their permutation statistics are 2 x B x rows float64


def _cpdb(L, R):
    """CellPhoneDB's score as LIANA+ computes it for the observed labels."""
    s = (L + R) / 2
    return np.where((L == 0) | (R == 0), 0.0, s)


def _cellchat(L, R):
    from liana.method.sc._cellchat import cellchat
    p = L * R
    return p / (cellchat._kh + p)


def _connectome(L, R):
    """Connectome's magnitude as LIANA+ computes it (``expr_prod``)."""
    return L * R


SCORES = {"CellPhoneDB": (_cpdb, ("ligand_means", "receptor_means")),
          "CellChat": (_cellchat, ("ligand_trimean", "receptor_trimean")),
          "Connectome": (_connectome, ("ligand_means", "receptor_means"))}


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
        # LIANA+'s scores and ours are row by row, so rows are scored in blocks (see lazy_permutations)
        lazy = isinstance(perm_stats, _PermStats)
        blocks = [slice(a, a + ROWS_AT_ONCE) for a in range(0, len(x), ROWS_AT_ONCE)] if lazy else [slice(None)]
        out = []
        for z in blocks:
            xz, ps = x.iloc[z], perm_stats.rows(z) if lazy else perm_stats
            mag, p = method.fun(xz, ps)
            nb = side_p(score, xz[lcol].to_numpy(), xz[rcol].to_numpy(), ps)
            out.append((np.asarray(mag), np.asarray(p), nb))
        mag, p = (np.concatenate([o[i] for o in out]) for i in (0, 1))
        nb = {k: np.concatenate([o[2][k] for o in out]) for k in out[0][2]}
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
    """LIANA+'s Connectome with permutations: its output (``expr_prod``, ``scaled_weight``) is LIANA+'s own,
    unchanged; ``store['Connectome']`` gets its rows, a joint-permutation p (LIANA+'s formula,
    #{s_b >= s} / B, on ``expr_prod``) and Null B on the same permutations."""
    from liana.method._pipe_utils._get_mean_perms import _calculate_pvals
    from liana.method.sc._Method import Method, MethodMeta

    method = __import__("liana").mt.connectome
    score, (lcol, rcol) = SCORES["Connectome"]

    def fun(x, perm_stats):
        mag, spec = method.fun(x)
        lazy = isinstance(perm_stats, _PermStats)
        blocks = [slice(a, a + ROWS_AT_ONCE) for a in range(0, len(x), ROWS_AT_ONCE)] if lazy else [slice(None)]
        pa, nbs = [], []
        for z in blocks:
            xz, ps = x.iloc[z], np.asarray(perm_stats.rows(z) if lazy else perm_stats)
            pa.append(_calculate_pvals(np.asarray(mag)[z], ps, np.prod))
            nbs.append(side_p(score, xz[lcol].to_numpy(), xz[rcol].to_numpy(), ps))
        nb = {k: np.concatenate([o[k] for o in nbs]) for k in nbs[0]}
        store["Connectome"] = pd.DataFrame({**{k: x[k].to_numpy() for k in KEY + ["ligand", "receptor"]},
                                            "magnitude": np.asarray(mag), "p_null_a": np.concatenate(pa),
                                            "p_L_null_b": nb["p_L"], "p_R_null_b": nb["p_R"], "p_null_b": nb["p"],
                                            "score_check": nb["score"]})
        return mag, spec

    meta = MethodMeta(method_name=method.method_name, complex_cols=method.complex_cols, add_cols=method.add_cols,
                      fun=fun, magnitude=method.magnitude, magnitude_ascending=method.magnitude_ascending,
                      specificity=method.specificity, specificity_ascending=method.specificity_ascending,
                      permute=True, reference=method.reference)
    return Method(_method=meta)


class _PermStats:
    """LIANA+'s ``perm_stats`` (2, B, rows) before it is built: ``rows(z)`` builds rows ``z`` only."""

    def __init__(self, cube, sides):
        self.cube, self.sides = cube, sides          # sides: [(cell types, genes)] for ligand, receptor
        self.shape = (2, cube.shape[0], len(sides[0][0]))

    def rows(self, z):
        return np.stack([self.cube[:, c[z], g[z]] for c, g in self.sides], axis=0)

    def __array__(self, dtype=None, copy=None):      # any other consumer gets the whole array, as in LIANA+
        a = self.rows(slice(None))
        return a if dtype is None else a.astype(dtype)


class _Side:
    """``perms[:, cell types, genes]`` not yet taken; ``np.stack`` of the two sides gives a ``_PermStats``."""

    def __init__(self, cube, c, g):
        self.cube, self.c, self.g = cube, c, g

    def __array_function__(self, func, types, args, kwargs):
        if func is np.stack and kwargs.get("axis", 0) == 0 and all(isinstance(a, _Side) for a in args[0]):
            return _PermStats(self.cube, [(a.c, a.g) for a in args[0]])
        return func(*[np.asarray(self.cube[:, self.c, self.g]) if a is self else a for a in args], **kwargs)

    def __array__(self, dtype=None, copy=None):
        a = self.cube[:, self.c, self.g]
        return a if dtype is None else a.astype(dtype)


class _Cube:
    def __init__(self, a):
        self.a, self.shape = a, a.shape

    def __getitem__(self, key):
        if isinstance(key, tuple) and len(key) == 3 and key[0] == slice(None):
            return _Side(self.a, np.asarray(key[1]), np.asarray(key[2]))
        return self.a[key]


@contextmanager
def lazy_permutations():
    """Within this block LIANA+ does not build ``perm_stats`` for every row at once: the captured tools
    score rows in blocks of ``ROWS_AT_ONCE``. Same permutations, same arithmetic per row, so the same
    numbers; memory is 2 x B x ROWS_AT_ONCE floats instead of 2 x B x rows (over 30 GB on Tabula
    Sapiens Lung: 34 cell types)."""
    import liana.method.sc._liana_pipe as pipe

    inner = pipe._get_means_perms
    pipe._get_means_perms = lambda *a, **k: _Cube(inner(*a, **k))
    try:
        yield
    finally:
        pipe._get_means_perms = inner


@contextmanager
def stratified_permutations(strata: pd.Series):
    """Within this block, LIANA+ permutes cell labels only within each stratum (``strata``: one value
    per cell, indexed by obs name, e.g. the donor). Same seed handling and aggregation as LIANA+'s
    ``_get_means_perms``; every cell must have a stratum."""
    import liana.method.sc._liana_pipe as pipe
    from liana.method._pipe_utils import _get_mean_perms as gm

    original = pipe._get_means_perms

    def perms(adata, n_perms, seed, agg_fun, norm_factor, n_jobs, verbose):
        if isinstance(norm_factor, np.float32):
            adata.X /= norm_factor
        labels = adata.obs["@label"].cat.categories
        mask = np.stack([(adata.obs["@label"] == lab).to_numpy() for lab in labels], axis=1)
        s = strata.reindex(adata.obs_names)
        if s.isna().any():
            raise ValueError(f"{int(s.isna().sum())} cells have no stratum")
        codes = pd.Categorical(s.astype(str)).codes
        groups = [np.where(codes == k)[0] for k in np.unique(codes)]
        rng = np.random.default_rng(seed)
        out = np.zeros((n_perms, len(labels), adata.n_vars))
        for b in range(n_perms):
            idx = np.arange(adata.n_obs)
            for g in groups:
                idx[g] = g[rng.permutation(len(g))]
            _, m = gm._permute_and_aggregate(b, idx, adata.X, mask, agg_fun)
            out[b] = np.reshape(m, (len(labels), adata.n_vars))
        return out

    pipe._get_means_perms = perms
    try:
        yield
    finally:
        pipe._get_means_perms = original


def run(adata, groupby: str, *, resource_name: str = "consensus", n_perms: int = 1000, seed: int = 1337,
        expr_prop: float = 0.1, strata: pd.Series | None = None, tools_only: bool = False, use_raw: bool = False,
        connectome: bool = False, **kw) -> tuple[dict, dict]:
    """One LIANA+ run of the six methods with CellPhoneDB and CellChat captured, and Connectome too if
    ``connectome`` (``tools_only``: just the captured ones). Returns ({method name: its table, keys +
    magnitude + specificity}, {tool: captured rows})."""
    import liana as li
    from liana.method.sc._rank_aggregate import AggregateClass, _rank_aggregate_meta

    store = {}
    methods = [capture(li.mt.cellphonedb, store), capture_connectome(store) if connectome else li.mt.connectome,
               li.mt.logfc, li.mt.natmi, li.mt.singlecellsignalr,
               capture(li.mt.cellchat, store)]   # CellChat last: it rescales X in place
    if tools_only:
        methods = [methods[0]] + ([methods[1]] if connectome else []) + [methods[-1]]
    agg = AggregateClass(_rank_aggregate_meta, methods=methods)
    kw = dict(groupby=groupby, resource_name=resource_name, n_perms=n_perms, seed=seed, expr_prop=expr_prop,
              consensus_opts=False, inplace=False, use_raw=use_raw, **kw)
    if strata is None:
        with lazy_permutations():
            lrs = agg(adata, **kw)
    else:
        with stratified_permutations(strata), lazy_permutations():
            lrs = agg(adata, **kw)
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


def variants(lrs: dict, store: dict, nulls=NULLS, aggregate_method: str = "rra") -> dict:
    """Every consensus variant, each a LIANA+ consensus table (keys, the methods' scores,
    ``magnitude_rank``, ``specificity_rank``)."""
    import liana as li
    from liana.method._pipe_utils._aggregate import _aggregate
    from liana.method.sc._rank_aggregate import AggregateClass, _rank_aggregate_meta

    six = {m.method_name: m for m in [li.mt.cellphonedb, li.mt.connectome, li.mt.logfc, li.mt.natmi,
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
    if "Connectome" in store:
        metas3 = [_meta("Connectome", *CONNECTOME) if m.method_name == "Connectome" else m for m in six.values()]
        for n in nulls:
            swapped = dict(base)
            for tool in TOOLS if n != "null_a" else ():
                swapped[tool] = ours(tool, n)
            d = store["Connectome"][KEY + ["magnitude", f"p_{n}"]].drop_duplicates(KEY)
            swapped["Connectome"] = d.rename(columns={"magnitude": CONNECTOME[0], f"p_{n}": CONNECTOME[1]})
            out[f"swap3_{n}"] = agg(swapped, metas3)
    return out


DONOR_KEYS = ("donor_id", "donor", "patient_id", "patient", "individual", "subject", "orig.ident", "sample_id", "sample")


def donor_strata(adata, groupby: str, strata: str | None = "auto") -> str | None:
    """The ``obs`` column whose values the permutations stay within, or None for LIANA+'s global permutations.

    ``"auto"``: the first of ``DONOR_KEYS`` in ``adata.obs`` with two or more values, provided at least two of its
    groups hold two or more cell types (else no label could move and the test would have no power: global
    permutations, with a warning). A column name is used as given; None means global permutations.

    Why: with several donors and labels that follow donors, global permutations let donor differences pass as
    cell-type specificity. On labels shuffled within donor, Null B with global permutations gave p <= 0.05 for
    11-18% of rows (TNBC, Tabula Sapiens Liver), and 0.3% with permutations within donor (AMENDMENT_2, A2.3).
    """
    import warnings
    if strata != "auto":
        if strata is not None and strata not in adata.obs:
            raise KeyError(f"strata: no column {strata!r} in adata.obs")
        return strata
    for k in DONOR_KEYS:
        if k in adata.obs and adata.obs[k].nunique() >= 2:
            types_per_group = adata.obs.groupby(adata.obs[k].astype(str), observed=True)[groupby].nunique()
            if (types_per_group >= 2).sum() < 2:
                warnings.warn(f"strata='auto': {k!r} has {adata.obs[k].nunique()} values but fewer than two hold two "
                              f"or more cell types; using LIANA+'s global permutations", stacklevel=3)
                return None
            return k
    return None


def rank_aggregate(adata, groupby: str, *, null: str = "null_b", mode: str = "swap", resource_name: str = "consensus",
                   n_perms: int = 1000, seed: int = 1337, expr_prop: float = 0.1, strata: str | None = "auto",
                   aggregate_method: str = "rra", use_raw: bool = False, key_added: str = "liana_res",
                   inplace: bool = True, return_tools: bool = False, **kw):
    """LIANA+'s ``rank_aggregate`` (six methods: its five plus CellChat) with the permutation tools'
    p-values from our null.

    Parameters
    ----------
    adata
        Log-normalised expression in ``X``; it is not modified.
    groupby
        Column of ``adata.obs`` with the cell types.
    null
        ``"null_b"`` (composite null: both sides must be specific) or ``"null_a"`` (LIANA+'s own p).
    mode
        ``"swap"``: CellPhoneDB's and CellChat's p replaced by the null's (recommended);
        ``"extra"``: the null's p added as two more methods (7th and 8th).
    strata
        ``"auto"`` (default): with several donors, permutations stay within each donor. The donor column is the
        first of ``DONOR_KEYS`` present with two or more values (see :func:`donor_strata`); without one, LIANA+'s
        global permutations. A column name: permutations stay within each of its values. None: global
        permutations, as LIANA+. The column used is in the result's ``attrs["strata"]``.
    aggregate_method
        ``"rra"`` (LIANA+'s default) or ``"mean"``.
    use_raw
        Score ``adata.raw`` instead of ``X``.
    return_tools
        Also return the two tools' rows with every p (``p_null_a``, ``p_L_null_b``, ``p_R_null_b``, ``p_null_b``).

    Returns
    -------
    The consensus table (``magnitude_rank``, ``specificity_rank`` and each method's columns), also written to
    ``adata.uns[key_added]`` when ``inplace``; with ``return_tools``, ``(table, tools)``.

    Examples
    --------
    >>> import scanpy as sc
    >>> from ccc_nulls import liana_consensus as lc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> res = lc.rank_aggregate(adata, "bulk_labels", use_raw=True, n_perms=100, inplace=False)
    """
    if null not in ("null_a", "null_b"):
        raise ValueError("null must be 'null_a' or 'null_b' (B2 needs decoy genes: see add_b2)")
    if mode not in ("swap", "extra"):
        raise ValueError("mode must be 'swap' or 'extra'")
    strata = donor_strata(adata, groupby, strata)
    s = adata.obs[strata].astype(str) if strata is not None else None
    lrs, store = run(adata.copy(), groupby, resource_name=resource_name, n_perms=n_perms, seed=seed,
                     expr_prop=expr_prop, strata=s, use_raw=use_raw, **kw)
    V = variants(lrs, store, nulls=("null_a", "null_b"), aggregate_method=aggregate_method)
    name = "liana6" if (null, mode) == ("null_a", "swap") else f"{mode}_{null}"
    res = V[name].sort_values("specificity_rank", kind="stable").reset_index(drop=True)
    res.attrs["strata"] = strata
    if inplace:
        adata.uns[key_added] = res
    tools = {k: d.drop(columns="score_check") for k, d in store.items()}
    return (res, tools) if return_tools else res
