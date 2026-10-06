"""Thin adapters from native tool output to the common schema.

Common schema: ``ligand, receptor, source, target, tool_score`` (one row per gene-level
candidate), plus optional ``ligand_complex`` / ``receptor_complex`` for
:func:`ccc_nulls.core.aggregate_complexes`.

Native CellPhoneDB / CellChat (R) adapters are on the roadmap. Until then, run those
methods through liana-py and use :func:`from_liana`, or use :func:`to_schema`.
"""
from __future__ import annotations

import pandas as pd

# liana-py magnitude column per method (checked against liana 1.x: li.mt.<method>.magnitude;
# logfc has no magnitude column, so its specificity score is used)
LIANA_SCORES = {
    "cellphonedb": "lr_means",
    "connectome": "expr_prod",
    "singlecellsignalr": "lrscore",
    "natmi": "expr_prod",
    "logfc": "lr_logfc",
    "cellchat": "lr_probs",
    "geometric_mean": "lr_gmeans",
    "scseqcomm": "inter_score",
    "rank_aggregate": "magnitude_rank",
}


def to_schema(df: pd.DataFrame, ligand: str, receptor: str, source: str, target: str,
              tool_score: str | None = None, keep: list[str] | None = None) -> pd.DataFrame:
    """Generic adapter: rename columns into the common schema."""
    cols = {ligand: "ligand", receptor: "receptor", source: "source", target: "target"}
    if tool_score:
        cols[tool_score] = "tool_score"
    out = df.rename(columns=cols)[list(cols.values()) + (keep or [])].copy()
    return out.reset_index(drop=True)


def from_liana(res: pd.DataFrame, method: str | None = None, tool_score: str | None = None) -> pd.DataFrame:
    """Convert a liana-py result table (``adata.uns[key_added]``).

    liana-py reports one row per (complex) interaction. For a heteromeric complex, the
    ``ligand`` / ``receptor`` column holds the representative subunit (the lowest-expressed
    one, as in CellPhoneDB's minimum-subunit rule), and that subunit is what both nulls
    test. ``ligand_complex`` / ``receptor_complex`` are kept for reporting.
    """
    if tool_score is None and method is not None:
        tool_score = LIANA_SCORES.get(method.lower())
    if tool_score is not None and tool_score not in res.columns:
        raise KeyError(f"column {tool_score!r} not in liana result; columns are {list(res.columns)}")
    keep = [c for c in ("ligand_complex", "receptor_complex") if c in res.columns]
    out = to_schema(res, "ligand", "receptor", "source", "target", tool_score, keep)
    return out.drop_duplicates(["ligand", "receptor", "source", "target"] + keep).reset_index(drop=True)

