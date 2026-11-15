"""Dense / scipy-sparse helpers. Everything downstream works on cells x genes matrices."""
from __future__ import annotations

import numpy as np
from scipy import sparse


def as_matrix(X):
    """Return X as a float CSC matrix (if sparse) or a float ndarray."""
    if sparse.issparse(X):
        return sparse.csc_matrix(X, dtype=np.float64)
    return np.asarray(X, dtype=np.float64)


def onehot(labels: np.ndarray, n_types: int) -> np.ndarray:
    """types x cells averaging matrix (each row sums to 1)."""
    O = np.zeros((n_types, len(labels)))
    O[labels, np.arange(len(labels))] = 1.0
    sizes = O.sum(1, keepdims=True)
    if np.any(sizes == 0):
        raise ValueError("every cell type must contain at least one cell")
    return O / sizes


def cluster_means(X, labels: np.ndarray, n_types: int) -> np.ndarray:
    """types x genes matrix of mean expression per cell type."""
    O = onehot(labels, n_types)
    if sparse.issparse(X):
        return np.asarray((X.T @ O.T).T)
    return O @ X


def gene_mean_var(X) -> tuple[np.ndarray, np.ndarray]:
    """Per-gene mean and population variance (ddof = 0)."""
    if sparse.issparse(X):
        mu = np.asarray(X.mean(0)).ravel()
        sq = np.asarray(X.multiply(X).mean(0)).ravel()
        return mu, np.maximum(sq - mu**2, 0.0)
    return X.mean(0), X.var(0)


def detection(X) -> np.ndarray:
    """Per-gene fraction of cells with non-zero expression."""
    if sparse.issparse(X):
        return np.asarray((X > 0).mean(0)).ravel()
    return (X > 0).mean(0)


def detection_by_type(X, labels: np.ndarray, n_types: int) -> np.ndarray:
    """types x genes fraction of cells expressing each gene (counts / size, so exact
    ratios such as 50/500 compare cleanly against a 0.10 threshold)."""
    ind = np.zeros((n_types, len(labels)))
    ind[labels, np.arange(len(labels))] = 1.0
    nz = (X > 0).astype(np.float64)
    counts = np.asarray((nz.T @ ind.T).T) if sparse.issparse(X) else ind @ nz
    return counts / ind.sum(1, keepdims=True)


def columns(X, idx: np.ndarray):
    """Column subset that keeps the sparse / dense type."""
    return X[:, idx]
