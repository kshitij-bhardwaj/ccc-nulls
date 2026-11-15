"""ccc-nulls: a composite null for ligand-receptor specificity inside LIANA+'s consensus, and calibrated
nulls for cell–cell communication inference.

>>> import ccc_nulls as cn
>>> res = cn.rank_aggregate(adata, "cell_type")                  # LIANA+'s consensus with Null B (needs liana)
>>> res = cn.significance_layer(X, genes, labels, candidates)   # adds p/q under Null A and Null B
>>> cn.overlap_summary(res)                                      # both / A_only / B_only / none
>>> cn.composite.null_b_table(X, genes, labels, candidates)     # Null B as a composite null (both sides)
>>> cn.competitive.null_b2(X, gs, codes, rows, excluded)         # B2: decoy genes, CAMERA-corrected
"""
from . import competitive, composite, liana_consensus, simulate
from .core import aggregate_complexes, flag_overlap, overlap_summary, significance_layer
from .liana_consensus import rank_aggregate
from .nulls import gene_strata, null_a, null_b, null_b_binned_labels, perm_moments
from .scores import SCORES, score_mean, score_product
from .stats import bh, cauchy_combine

__version__ = "0.2.0"

__all__ = [
    "rank_aggregate", "liana_consensus",
    "significance_layer",
    "overlap_summary", "flag_overlap", "aggregate_complexes",
    "null_a", "null_b", "null_b_binned_labels", "gene_strata", "perm_moments",
    "score_mean", "score_product", "SCORES", "bh", "cauchy_combine",
    "competitive", "composite", "simulate",
]
