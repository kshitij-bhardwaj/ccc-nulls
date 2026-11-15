"""Open Problems v1 cell-cell communication metrics, reproduced exactly.

Works for both tasks: ligand-target (keys ``["ligand", "target"]``, cytokine truth) and
source-target (keys ``["source", "target"]``, spatial truth). Source: openproblems v1.0.0,
tasks/_cell_cell_communication/_common/{metrics/auprc.py, metrics/odds_ratio.py, utils.py}.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def predictions(rows: pd.DataFrame, score, keys, p=None, alpha: float = 0.05, how: str = "max",
                ligand_col: str = "ligand_complex") -> pd.DataFrame:
    """One score per truth key.

    If ``p`` is given, a row keeps its score only when p <= alpha (else 0), which is how
    Open Problems uses CellPhoneDB's permutation p-value. Rows are then aggregated per key.
    """
    s = np.asarray(score, dtype=float)
    if p is not None:
        s = np.where(np.asarray(p) <= alpha, s, 0.0)
    cols = {k: rows[ligand_col if k == "ligand" else k].to_numpy() for k in keys}
    return pd.DataFrame({**cols, "score": s}).groupby(list(keys)).score.agg(how).reset_index()


def _join(pred: pd.DataFrame, truth: pd.DataFrame, keys) -> pd.DataFrame:
    gt = truth.merge(pred, on=list(keys), how="left")
    gt["response"] = gt["response"].fillna(0)
    if gt["score"].notna().any():
        gt["score"] = gt["score"].fillna(np.nanmin(gt["score"]) - np.finfo(float).eps)
    else:
        gt["score"] = 0.0
    return gt


def auprc(pred: pd.DataFrame, truth: pd.DataFrame, keys) -> float:
    from sklearn.metrics import auc, precision_recall_curve

    gt = _join(pred, truth, keys)
    precision, recall, _ = precision_recall_curve(gt["response"], gt["score"], pos_label=1)
    return float(auc(recall, precision))


def odds_ratio(pred: pd.DataFrame, truth: pd.DataFrame, keys, top_prop: float = 0.05) -> float:
    """Raw odds ratio of positives in the top ``top_prop`` of truth rows vs the rest."""
    gt = _join(pred, truth, keys).sort_values("score", ascending=False)
    top = np.zeros(len(gt), bool)
    top[: int(truth.shape[0] * top_prop)] = True
    r = gt["response"].to_numpy() == 1
    tp, fp, fn, tn = (top & r).sum(), (top & ~r).sum(), (~top & r).sum(), (~top & ~r).sum()
    if fp * fn == 0:
        return float("nan") if tp * tn == 0 else float("inf")
    return float(tp * tn / (fp * fn))


def metrics(pred: pd.DataFrame, truth: pd.DataFrame, keys, top_prop: float = 0.05) -> dict:
    """AUPRC, odds ratio, and Open Problems' sigmoid-scaled odds ratio."""
    orr = odds_ratio(pred, truth, keys, top_prop)
    sig = 1 - 1 / (1 + orr / 2) if np.isfinite(orr) else (1.0 if orr == np.inf else float("nan"))
    return dict(auprc=auprc(pred, truth, keys), odds_ratio=orr, odds_ratio_sigmoid=float(sig))


def random_baseline(truth: pd.DataFrame, keys, seed: int = 1) -> dict:
    """Open Problems' "Random Events" control."""
    rnd = truth[list(keys)].copy()
    rnd["score"] = np.random.default_rng(seed).random(len(rnd))
    return metrics(rnd, truth, keys)
