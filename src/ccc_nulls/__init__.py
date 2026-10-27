"""ccc-nulls: structure-preserving null models and FDR control for cell–cell communication inference.

>>> import ccc_nulls as cn
>>> res = cn.significance_layer(X, genes, labels, candidates)   # adds p/q under Null A and Null B
>>> cn.overlap_summary(res)                                      # both / A_only / B_only / none
>>> cn.composite.null_b_table(X, genes, labels, candidates)     # Null B as a composite null (both sides)
"""
from . import adapters, composite, simulate
from .core import aggregate_complexes, flag_overlap, overlap_summary, significance_layer
from .io import from_anndata, significance_layer_anndata
from .nulls import gene_strata, null_a, null_b, null_b_binned_labels, perm_moments
from .scores import SCORES, score_mean, score_product
from .stats import bh, cauchy_combine

__version__ = "0.1.0"

__all__ = [
    "significance_layer", "significance_layer_anndata", "from_anndata",
    "overlap_summary", "flag_overlap", "aggregate_complexes",
    "null_a", "null_b", "null_b_binned_labels", "gene_strata", "perm_moments",
    "score_mean", "score_product", "SCORES", "bh", "cauchy_combine",
    "adapters", "composite", "simulate",
]
