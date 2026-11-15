"""LIANA's re-implementation of each tool, returned in the common schema.

Every tool workspace uses this as its reference: a native adapter should agree with it
wherever LIANA mirrors the original tool, and differences are worth documenting.
"""
from __future__ import annotations

import numpy as np

from . import schema

# liana-py method -> (statistic tested, native p-value column or None)
METHODS = {
    "cellphonedb": ("lr_means", "cellphone_pvals"),
    "cellchat": ("lr_probs", "cellchat_pvals"),
    "natmi": ("spec_weight", None),
    "connectome": ("scaled_weight", None),
    "geometric_mean": ("lr_gmeans", "gmean_pvals"),
    "singlecellsignalr": ("lrscore", None),
    "logfc": ("lr_logfc", None),
}


def run(adata, method: str, groupby: str, *, expr_prop: float = 0.1, n_perms: int | None = 1000,
        resource_name: str = "consensus", use_raw: bool = False, layer: str | None = None, seed: int = 0):
    """Run ``liana.mt.<method>`` and return the common schema."""
    import liana as li

    score_col, p_col = METHODS[method]
    kw = dict(groupby=groupby, resource_name=resource_name, expr_prop=expr_prop, use_raw=use_raw,
              layer=layer, verbose=False, key_added=f"liana_{method}", seed=seed)
    if p_col is not None:
        kw["n_perms"] = n_perms
    getattr(li.mt, method)(adata, **kw)
    res = adata.uns.pop(f"liana_{method}")
    out = res[["ligand", "receptor", "source", "target", "ligand_complex", "receptor_complex"]].copy()
    out["score"] = res[score_col].to_numpy(dtype=float)
    out["score_name"] = score_col
    out["p_native"] = res[p_col].to_numpy(dtype=float) if (p_col and p_col in res) else np.nan
    out["tool"] = method
    for c in ("source", "target"):
        out[c] = out[c].astype(str)
    return schema.validate(out)
