"""AnnData convenience wrappers (anndata is optional and imported lazily)."""
from __future__ import annotations

import pandas as pd

from .core import significance_layer


def from_anndata(adata, groupby: str, layer: str | None = None, use_raw: bool = False):
    """Return ``(X, genes, labels)`` from an AnnData object.

    Use the matrix your CCC tool scored. For liana-py with ``use_raw=True``, that is
    ``adata.raw``. The matrix must be normalised (log1p CP10k or similar), not scaled.
    """
    if use_raw:
        if adata.raw is None:
            raise ValueError("use_raw=True but adata.raw is None")
        X, genes = adata.raw.X, adata.raw.var_names
    else:
        X = adata.layers[layer] if layer else adata.X
        genes = adata.var_names
    return X, pd.Index(genes).astype(str), adata.obs[groupby].astype(str).to_numpy()


def significance_layer_anndata(adata, groupby: str, candidates: pd.DataFrame, *,
                               layer: str | None = None, use_raw: bool = False, **kwargs) -> pd.DataFrame:
    """:func:`ccc_nulls.significance_layer` taking an AnnData object."""
    X, genes, labels = from_anndata(adata, groupby, layer, use_raw)
    return significance_layer(X, genes, labels, candidates, **kwargs)
