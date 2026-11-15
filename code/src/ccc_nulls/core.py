"""Public entry point: add calibrated p/q-values under both nulls to any tool's output."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from ._matrix import as_matrix, detection_by_type
from .nulls import gene_strata, null_a, null_b
from .stats import bh

SCHEMA = ["ligand", "receptor", "source", "target"]


def significance_layer(
    X,
    genes,
    labels,
    candidates: pd.DataFrame,
    *,
    score="mean",
    B: int = 1000,
    alpha: float = 0.05,
    min_frac: float = 0.10,
    resource_genes=None,
    n_mean_bins: int = 10,
    n_det_bins: int = 3,
    max_decoys: int = 400,
    seed: int = 0,
) -> pd.DataFrame:
    """Recompute significance for a ranked LR table under Null A and Null B.

    Parameters
    ----------
    X
        cells x genes expression matrix: the *same* normalised (not scaled) matrix the
        tool scored. Dense ndarray or scipy sparse.
    genes, labels
        Gene names (length n_genes) and cell-type label per cell (length n_cells).
    candidates
        Table in the common schema, with columns ``ligand, receptor, source, target``
        plus anything else (kept as-is). Use :mod:`ccc_nulls.liana_reference` to convert tool
        output.
    score
        ``"mean"`` (CellPhoneDB), ``"product"`` (SingleCellSignalR-like), or a callable
        ``f(M, li, ri, si, ti)``. The score is recomputed from X, so both nulls test the
        same statistic.
    B
        Label permutations for Null A. Null B is enumerated exactly and needs no B.
    alpha
        q-value cutoff for the ``flag`` column.
    min_frac
        Expression filter the candidates passed. It is applied to Null B's decoys too.
    resource_genes
        All genes in the LR resource. They are excluded from the decoy pool. Default:
        the candidate genes.

    Returns
    -------
    DataFrame: the input rows that map to X, plus ``score_obs, p_rand, q_rand,
    z_strat, p_strat, q_strat, p_strat_min`` (the smallest p Null B can reach for that
    row, which is its resolution) and ``flag``, which is one of ``"both"``, ``"A_only"``,
    ``"B_only"`` or ``"none"``. ``A_only`` marks calls that are significant only
    because the uniform null ignored expression abundance.
    """
    missing = [c for c in SCHEMA if c not in candidates.columns]
    if missing:
        raise ValueError(f"candidates is missing columns {missing}; expected {SCHEMA}")
    X = as_matrix(X)
    genes = pd.Index(np.asarray(genes).astype(str))
    if X.shape[1] != len(genes):
        raise ValueError(f"X has {X.shape[1]} columns but {len(genes)} gene names were given")
    labels = np.asarray(labels).astype(str)
    if X.shape[0] != len(labels):
        raise ValueError(f"X has {X.shape[0]} rows but {len(labels)} labels were given")
    types = pd.Index(sorted(set(labels)))
    lab = types.get_indexer(labels)

    c = candidates.copy()
    li = genes.get_indexer(c.ligand.astype(str))
    ri = genes.get_indexer(c.receptor.astype(str))
    si = types.get_indexer(c.source.astype(str))
    ti = types.get_indexer(c.target.astype(str))
    ok = (li >= 0) & (ri >= 0) & (si >= 0) & (ti >= 0)
    if not ok.all():
        warnings.warn(f"dropping {int((~ok).sum())} candidates whose genes or cell types are not in X", stacklevel=2)
    c, li, ri, si, ti = c.loc[ok].copy(), li[ok], ri[ok], si[ok], ti[ok]
    if len(c) == 0:
        raise ValueError("no candidate maps onto X")

    decoy_pool = np.ones(len(genes), bool)
    decoy_pool[np.concatenate([li, ri])] = False
    if resource_genes is not None:
        rg = genes.get_indexer(pd.Index(np.asarray(list(resource_genes)).astype(str)))
        decoy_pool[rg[rg >= 0]] = False

    frac = detection_by_type(X, lab, len(types))
    low = (frac[si, li] < min_frac) | (frac[ti, ri] < min_frac)
    if low.any():
        warnings.warn(f"{int(low.sum())} candidates fall below min_frac={min_frac}; their decoys "
                      "are still filtered at min_frac, which makes them conservative", stacklevel=2)

    strata = gene_strata(X, n_mean_bins, n_det_bins)
    s_obs, p_rand = null_a(X, lab, len(types), li, ri, si, ti, score, B, seed)
    z, p_strat, n_ref = null_b(X, lab, len(types), li, ri, si, ti, score, seed + 1, strata=strata,
                               decoy_pool=decoy_pool, min_frac=min_frac, max_decoys=max_decoys,
                               return_n_ref=True)
    c["score_obs"], c["p_rand"], c["z_strat"], c["p_strat"] = s_obs, p_rand, z, p_strat
    c["p_strat_min"] = 1.0 / n_ref
    floor = alpha / len(c)                     # BH needs p <= alpha * k / m; k = 1 is the hardest
    if np.median(c.p_strat_min) > alpha * 50 / len(c):
        warnings.warn(f"Null B resolution is coarse: median smallest attainable p is "
                      f"{np.median(c.p_strat_min):.1e} against a BH threshold of {floor:.1e} "
                      f"for the top hit. Use more genes (the full transcriptome, not only HVGs) or "
                      f"fewer strata (n_mean_bins / n_det_bins)", stacklevel=2)
    c["q_rand"], c["q_strat"] = bh(p_rand), bh(p_strat)
    c["flag"] = flag_overlap(c.q_rand.to_numpy(), c.q_strat.to_numpy(), alpha)
    return c


def flag_overlap(q_rand, q_strat, alpha=0.05) -> np.ndarray:
    """Label each candidate by which null calls it significant at ``alpha``."""
    a, b = np.asarray(q_rand) < alpha, np.asarray(q_strat) < alpha
    return np.select([a & b, a & ~b, ~a & b], ["both", "A_only", "B_only"], "none")


def overlap_summary(res: pd.DataFrame) -> pd.Series:
    """Counts of each flag, the headline of every Null A vs Null B comparison."""
    return res["flag"].value_counts().reindex(["both", "A_only", "B_only", "none"], fill_value=0)


def aggregate_complexes(res: pd.DataFrame, ligand_col="ligand_complex", receptor_col="receptor_complex",
                        alpha: float = 0.05) -> pd.DataFrame:
    """Collapse subunit-level rows into complex-level rows, for tools whose adapters emit
    one row per subunit (liana-py output is already one row per complex).

    A complex is significant only if every subunit is (intersection-union test), so the
    complex p-value is the maximum over subunit rows. q-values and flags are
    recomputed at the complex level.
    """
    keys = [ligand_col, receptor_col, "source", "target"]
    agg = res.groupby(keys, sort=False).agg(p_rand=("p_rand", "max"), p_strat=("p_strat", "max"),
                                            score_obs=("score_obs", "min")).reset_index()
    agg["q_rand"], agg["q_strat"] = bh(agg.p_rand), bh(agg.p_strat)
    agg["flag"] = flag_overlap(agg.q_rand, agg.q_strat, alpha)
    return agg
