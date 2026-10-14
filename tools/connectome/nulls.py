"""Connectome: Null A and Null B on this tool's own score (``scaled_weight``).

Connectome has no permutation null on its edge score, so Null A here is our own: shuffle
labels and recompute the scaled edge weight (z-scores use the same scaling as the tool).

Both functions take the interactions the tool reported (common schema, from
``adapter.liana_reference`` or ``adapter.native_to_schema``) and the expression the tool scored,
and return the same rows, in the same order, with one p-value column added. The output is
checked by ``ccc_nulls.schema.validate_null``. The contract is in tools/README.md.

Shared helpers you may use (optional; ask before changing them, they live in src/ccc_nulls):
    ccc_nulls.multitool.StatsEngine     per-cell-type means / fractions / z-scores / trimeans
                                        under any labelling
    ccc_nulls.multitool.permute_codes   uniform or within-group (patient) label permutation
    ccc_nulls.nulls.gene_strata         mean-decile x detection-tercile gene strata
    ccc_nulls.stats.bh                  Benjamini-Hochberg q-values
"""
from __future__ import annotations

import pandas as pd


def null_a(adata, interactions: pd.DataFrame, groupby: str = "label", *, B: int = 1000,
           seed: int = 0) -> pd.DataFrame:
    """Null A: shuffle the cell-type labels, recompute Connectome's score, rank the observed score.

    Returns ``interactions`` with ``p_null_a`` added.
    """
    raise NotImplementedError("connectome: null_a is not implemented yet (see tools/connectome/nulls.py)")


def null_b(adata, interactions: pd.DataFrame, groupby: str = "label", *, B: int = 1000,
           seed: int = 0) -> pd.DataFrame:
    """Null B: keep the labels and rank the observed score against abundance-matched decoy genes.

    Returns ``interactions`` with ``p_null_b`` added.
    """
    raise NotImplementedError("connectome: null_b is not implemented yet (see tools/connectome/nulls.py)")
