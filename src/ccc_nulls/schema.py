"""The common interaction schema every tool adapter must produce.

One row per candidate interaction (ligand L in source cell type -> receptor R in target type):

    ligand, receptor          gene symbols actually tested (for a complex: the subunit used)
    source, target            cell-type labels, exactly as in the dataset's groupby column
    score                     the tool's own statistic (higher = stronger), float
    score_name                name of that statistic in the tool, e.g. "lr_means"
    p_native                  the tool's own p-value if it has one, else NaN
    ligand_complex,           complex identifiers (equal to ligand / receptor for single genes)
    receptor_complex
    tool                      tool identifier, e.g. "cellphonedb"

The null models in ``ccc_nulls`` consume this table plus the expression matrix the tool scored.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

KEY = ["ligand", "receptor", "source", "target"]
REQUIRED = KEY + ["score"]
OPTIONAL = ["score_name", "p_native", "ligand_complex", "receptor_complex", "tool"]
COLUMNS = REQUIRED + OPTIONAL


def validate(df: pd.DataFrame, *, require_p_native: bool = False) -> pd.DataFrame:
    """Return ``df`` in canonical column order, or raise ValueError listing every problem."""
    problems = []
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        problems.append(f"missing required columns {missing}")
    else:
        for c in KEY:
            if df[c].isna().any():
                problems.append(f"{int(df[c].isna().sum())} missing values in '{c}'")
        if not pd.api.types.is_numeric_dtype(df["score"]):
            problems.append("'score' must be numeric")
        elif not np.isfinite(df["score"].to_numpy(dtype=float)).all():
            problems.append("'score' has NaN or infinite values")
        keys = KEY + [c for c in ("ligand_complex", "receptor_complex") if c in df.columns]
        dup = int(df.duplicated(keys).sum())
        if dup:
            problems.append(f"{dup} duplicated rows on {keys}")
    if "p_native" in df.columns:
        p = df["p_native"].dropna()
        if len(p) and ((p < 0) | (p > 1)).any():
            problems.append("'p_native' outside [0, 1]")
    elif require_p_native:
        problems.append("'p_native' required for this tool")
    if problems:
        raise ValueError("schema check failed: " + "; ".join(problems))
    out = df.copy()
    for c in ("ligand_complex", "receptor_complex"):
        if c not in out.columns:
            out[c] = out[c.split("_")[0]]
    if "p_native" not in out.columns:
        out["p_native"] = np.nan
    extra = [c for c in out.columns if c not in COLUMNS]
    return out[[c for c in COLUMNS if c in out.columns] + extra].reset_index(drop=True)
