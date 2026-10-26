"""Simulator with known ground truth, for calibrating null models.

Regimes
-------
S0  ``structured=False``, labels shuffled afterwards: no cell-type structure, no CCC.
S1  ``structured=True, planted=False``: every gene carries a cell-type effect
    ``delta ~ N(0, tau^2)`` (log scale), but no LR pair is special. Every call is a
    false positive.
S2  ``structured=True, planted=True``: S1 plus ``n_planted`` events with the ligand up
    in X and the receptor up in Y by ``plant_lfc``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class SimData:
    X: np.ndarray                     # cells x genes, log1p(CP10k)
    genes: np.ndarray
    labels: np.ndarray                # cell-type names
    pairs: pd.DataFrame               # the LR "resource": ligand, receptor
    truth: set = field(default_factory=set)      # {(ligand, receptor, source, target)}
    lig_up: set = field(default_factory=set)     # {(ligand, source)} planted
    rec_up: set = field(default_factory=set)     # {(receptor, target)} planted


def simulate(seed=0, n_genes=8000, sizes=(500, 400, 300, 200, 150, 100), tau=0.25,
             n_pairs=400, n_planted=60, plant_lfc=1.0, dispersion=0.3,
             structured=True, planted=True, shuffle_labels=False) -> SimData:
    """Negative-binomial counts with a broad abundance spectrum, CP10k + log1p normalised."""
    rng = np.random.default_rng(seed)
    K, n = len(sizes), sum(sizes)
    codes = np.repeat(np.arange(K), sizes)
    log_mu = rng.normal(-1.0, 1.5, n_genes)
    delta = rng.normal(0, tau, (K, n_genes)) if structured else np.zeros((K, n_genes))
    pairs = np.array([rng.choice(n_genes, 2, replace=False) for _ in range(n_pairs)])
    truth, lig_up, rec_up = set(), set(), set()
    genes = np.array([f"G{i:05d}" for i in range(n_genes)])
    types = np.array([f"T{k}" for k in range(K)])
    if planted:
        for j in rng.choice(n_pairs, n_planted, replace=False):
            x, y = rng.choice(K, 2, replace=True)
            L, R = pairs[j]
            delta[x, L] += plant_lfc
            delta[y, R] += plant_lfc
            truth.add((genes[L], genes[R], types[x], types[y]))
            lig_up.add((genes[L], types[x]))
            rec_up.add((genes[R], types[y]))
    size_f = rng.lognormal(0, 0.3, n)
    mu = np.exp(log_mu[None, :] + delta[codes]) * size_f[:, None]
    counts = rng.poisson(rng.gamma(1 / dispersion, dispersion * mu))
    lib = np.maximum(counts.sum(1, keepdims=True), 1)
    X = np.log1p(counts / lib * 1e4)
    if shuffle_labels:
        codes = np.random.default_rng(seed + 99).permutation(codes)
    return SimData(X, genes, types[codes], pd.DataFrame({"ligand": genes[pairs[:, 0]], "receptor": genes[pairs[:, 1]]}),
                   truth, lig_up, rec_up)


def candidates_from_pairs(X, genes, labels, pairs: pd.DataFrame, min_frac=0.10) -> pd.DataFrame:
    """Every (L, R, X, Y) passing a CellPhoneDB-style expression filter, in schema form."""
    types = np.array(sorted(set(labels)))
    gi = pd.Index(genes)
    frac = np.stack([(X[labels == t] > 0).mean(0) for t in types])
    rows = []
    for L, R in zip(gi.get_indexer(pairs.ligand), gi.get_indexer(pairs.receptor)):
        for x in np.where(frac[:, L] >= min_frac)[0]:
            for y in np.where(frac[:, R] >= min_frac)[0]:
                rows.append((genes[L], genes[R], types[x], types[y]))
    return pd.DataFrame(rows, columns=["ligand", "receptor", "source", "target"])


def label_truth(cands: pd.DataFrame, sim: SimData) -> np.ndarray:
    """'true', 'partial' (one side planted), or 'null' for each candidate."""
    key = list(zip(cands.ligand, cands.receptor, cands.source, cands.target))
    true = np.array([k in sim.truth for k in key])
    part = np.array([(l, s) in sim.lig_up or (r, t) in sim.rec_up for l, r, s, t in key]) & ~true
    return np.where(true, "true", np.where(part, "partial", "null"))


# ------------------------------------------------------------------ named genes, real LR resource
COMPLEX_SEPARATOR = "_"        # LIANA writes heteromeric complexes as "SUB1_SUB2"


def load_pairs(path_or_df, *, simple_only: bool = True) -> pd.DataFrame:
    """Read an LR resource (e.g. LIANA's consensus CSV) as unique ``ligand, receptor`` pairs.

    Accepts ``ligand/receptor`` or ``source/target`` columns. With ``simple_only``,
    complexes (``"A_B"``) are dropped, so every planted event is a plain gene pair.
    """
    df = pd.read_csv(path_or_df) if not isinstance(path_or_df, pd.DataFrame) else path_or_df.copy()
    if {"ligand", "receptor"} <= set(df.columns):
        df = df[["ligand", "receptor"]]
    elif {"source", "target"} <= set(df.columns):
        df = df[["source", "target"]].set_axis(["ligand", "receptor"], axis=1)
    else:
        raise ValueError(f"resource needs ligand/receptor or source/target columns, got {list(df.columns)}")
    df = df.dropna().astype(str)
    if simple_only:
        df = df[~df.ligand.str.contains(COMPLEX_SEPARATOR) & ~df.receptor.str.contains(COMPLEX_SEPARATOR)]
    df = df[df.ligand != df.receptor].drop_duplicates().reset_index(drop=True)
    if df.empty:
        raise ValueError("resource has no usable ligand-receptor pairs")
    return df


def simulate_named(pairs: pd.DataFrame, seed: int = 0, n_background: int = 5000,
                   sizes=(800, 600, 500, 400, 400, 300), tau: float = 0.25, n_planted: int = 60,
                   plant_lfc: float = 1.0, dispersion: float = 0.3, structured: bool = True,
                   planted: bool = True, shuffle_labels: bool = False,
                   cell_types=None, background_prefix: str = "BG", plant_pairs=None) -> SimData:
    """The model of :func:`simulate`, but on a real LR resource so native tools can run it.

    LR genes carry their real symbols (from ``pairs``, e.g. LIANA's consensus via
    :func:`load_pairs`). ``n_background`` synthetic genes (``BG00000``...) fill out the
    transcriptome: no LR database knows them, so they act as decoys only. Every gene gets
    an abundance from the same distribution and, when ``structured``, a cell-type effect
    ``N(0, tau^2)``.

    Planted events use pairs whose ligand and receptor are not shared with any other
    planted event, so each planted gene marks exactly one true interaction. Other
    resource pairs that reuse a planted gene come out as ``"partial"`` in
    :func:`label_truth`.

    ``plant_pairs`` limits planting to pairs that are also in that table (same
    ligand -> receptor direction), e.g. the pairs every tool's own database contains.
    Then every planted event is visible to every tool. Background and resource genes are
    unchanged.

    Regimes as in :func:`simulate`: S0 ``structured=False, planted=False,
    shuffle_labels=True``. S1 ``planted=False``. S2 defaults.
    """
    pairs = load_pairs(pairs)
    pool = pairs if plant_pairs is None else pairs.merge(load_pairs(plant_pairs))
    if planted and pool.empty:
        raise ValueError("no resource pair is also in plant_pairs")
    rng = np.random.default_rng(seed)
    K, n = len(sizes), sum(sizes)
    types = np.array(cell_types if cell_types is not None else [f"CT{k}" for k in range(K)], dtype=str)
    if len(types) != K or len(set(types)) != K:
        raise ValueError("cell_types must give one unique name per entry of sizes")
    lr_genes = pd.unique(pairs[["ligand", "receptor"]].to_numpy().ravel())
    background = np.array([f"{background_prefix}{i:05d}" for i in range(n_background)])
    if set(background) & set(lr_genes):
        raise ValueError("background names clash with resource genes; change background_prefix")
    genes = np.concatenate([lr_genes, background])
    gi = pd.Index(genes)
    G = len(genes)

    codes = np.repeat(np.arange(K), sizes)
    log_mu = rng.normal(-1.0, 1.5, G)
    delta = rng.normal(0, tau, (K, G)) if structured else np.zeros((K, G))
    truth, lig_up, rec_up = set(), set(), set()
    if planted:
        used = set()
        for j in rng.permutation(len(pool)):
            if len(truth) == n_planted:
                break
            L, R = pool.ligand.iat[j], pool.receptor.iat[j]
            if L in used or R in used:
                continue
            used |= {L, R}
            x, y = rng.choice(K, 2, replace=True)
            delta[x, gi.get_loc(L)] += plant_lfc
            delta[y, gi.get_loc(R)] += plant_lfc
            x, y = str(types[x]), str(types[y])
            truth.add((L, R, x, y))
            lig_up.add((L, x))
            rec_up.add((R, y))
        if len(truth) < n_planted:
            raise ValueError(f"only {len(truth)} pairs with unshared genes for n_planted={n_planted}")
    size_f = rng.lognormal(0, 0.3, n)
    mu = np.exp(log_mu[None, :] + delta[codes]) * size_f[:, None]
    counts = rng.poisson(rng.gamma(1 / dispersion, dispersion * mu))
    lib = np.maximum(counts.sum(1, keepdims=True), 1)
    X = np.log1p(counts / lib * 1e4)
    if shuffle_labels:
        codes = np.random.default_rng(seed + 99).permutation(codes)
    return SimData(X, genes, types[codes], pairs, truth, lig_up, rec_up)


def restrict_truth(sim: SimData, resource: pd.DataFrame) -> set:
    """Planted events a given tool can see: those whose (ligand, receptor) is in its resource.

    Use it to score each tool against its own database. A planted pair missing from a
    tool's resource is not a miss for that tool.
    """
    known = set(map(tuple, load_pairs(resource)[["ligand", "receptor"]].to_numpy()))
    return {t for t in sim.truth if (t[0], t[1]) in known}


def sim_to_anndata(sim: SimData, **params):
    """SimData -> AnnData (sparse float32 X, ``obs["label"]``, truth tables in ``uns``)."""
    import anndata as ad
    from scipy import sparse

    def table(rows, cols):
        return pd.DataFrame(sorted(rows), columns=cols)

    a = ad.AnnData(X=sparse.csr_matrix(sim.X, dtype=np.float32),
                   obs=pd.DataFrame({"label": pd.Categorical(sim.labels)},
                                    index=[f"cell{i:05d}" for i in range(len(sim.labels))]),
                   var=pd.DataFrame(index=pd.Index(sim.genes, name=None)))
    a.uns["sim_truth"] = table(sim.truth, ["ligand", "receptor", "source", "target"])
    a.uns["sim_lig_up"] = table(sim.lig_up, ["ligand", "source"])
    a.uns["sim_rec_up"] = table(sim.rec_up, ["receptor", "target"])
    a.uns["sim_pairs"] = sim.pairs.reset_index(drop=True)
    a.uns["sim_params"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in params.items()}
    return a


def sim_from_anndata(adata) -> SimData:
    """Inverse of :func:`sim_to_anndata`, so :func:`label_truth` works on a saved simulation."""
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    u = adata.uns
    return SimData(np.asarray(X, dtype=float), np.asarray(adata.var_names).astype(str),
                   adata.obs["label"].astype(str).to_numpy(), pd.DataFrame(u["sim_pairs"]),
                   set(map(tuple, pd.DataFrame(u["sim_truth"]).to_numpy())),
                   set(map(tuple, pd.DataFrame(u["sim_lig_up"]).to_numpy())),
                   set(map(tuple, pd.DataFrame(u["sim_rec_up"]).to_numpy())))
