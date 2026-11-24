"""NicheNet (Browaeys et al., Nat Methods 2020; v2 matrix, Zenodo 7074291) as a downstream footprint.

NOT pre-registered: an addition requested after AMENDMENT_5's results were known.

Two uses, kept apart:

1. Truth (``weights`` + ``downstream.truth_rows(..., "nichenet")``): a third footprint resource next to CytoSig and
   PROGENy, scored exactly as they are (Amendment 1, A2): the receiver's contrast profile regressed on the
   ligand's footprint (decoupler's ulm), z-scored across receivers, responds if z >= 1.645. A ligand's footprint
   is its regulatory potential over its top ``TOP`` targets (the paper's cutoff for top targets), 0 elsewhere.

2. Filter p-value (``footprint_p``): for a ligand L and receiver t, is t's footprint of L stronger than its
   footprint of other ligands? Activities (ulm t-values) of every NicheNet ligand in every receiver are
   z-scored per ligand across receivers (a ligand whose targets are up everywhere does not stand out), then
   p(L, t) = rank of L among all ligands within t: (1 + #{L' != L : z(L', t) >= z(L, t)}) / n_ligands. A
   competitive null over ligands, not over genes, so correlated target genes (modules) do not shrink it.
   Resolution 1 / n_ligands (~1/1,000).

Data: data/nichenet/ligand_target_matrix.{bin,rows.txt,cols.txt}, converted from
ligand_target_matrix_nsga2r_final.rds (md5 b09606b04b2d4490418d9028c0e58b9f) by R's readRDS + writeBin (float32,
column-major).
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

from ccc_nulls import datasets

TOP = 250


def data_dir():
    return datasets.data_dir() / "nichenet"


@lru_cache(None)
def matrix() -> pd.DataFrame:
    """targets x ligands regulatory potential (human symbols)."""
    d = data_dir()
    rows = (d / "ligand_target_matrix.rows.txt").read_text().splitlines()
    cols = (d / "ligand_target_matrix.cols.txt").read_text().splitlines()
    v = np.fromfile(d / "ligand_target_matrix.bin", dtype=np.float32)
    return pd.DataFrame(v.reshape(len(cols), len(rows)).T, index=rows, columns=cols)


@lru_cache(None)
def weights(top: int = TOP) -> pd.DataFrame:
    """targets x ligands: regulatory potential over each ligand's ``top`` targets, 0 elsewhere."""
    M = matrix()
    A = M.to_numpy()
    out = np.zeros_like(A)
    idx = np.argpartition(-A, top, axis=0)[:top]
    cols = np.arange(A.shape[1])
    out[idx, cols] = A[idx, cols]
    return pd.DataFrame(out, index=M.index, columns=M.columns)


def ligand_of(ligand_complex: str, ligands) -> str | None:
    """The NicheNet ligand for a (human-symbol) ligand complex: its first subunit that NicheNet has."""
    return next((g for g in str(ligand_complex).split("_") if g in ligands), None)


def activities(profiles: pd.DataFrame, top: int = TOP) -> pd.DataFrame:
    """ligands x receivers: ulm slope t-value of each receiver's contrast profile on each ligand's footprint."""
    W = weights(top)
    genes = profiles.index.intersection(W.index)
    Y = profiles.loc[genes].to_numpy(float)
    X = W.loc[genes].to_numpy(float)
    n = len(genes)
    Xc, Yc = X - X.mean(0), Y - Y.mean(0)
    sxx = (Xc ** 2).sum(0)                         # ligands
    syy = (Yc ** 2).sum(0)                         # receivers
    b = (Xc.T @ Yc) / sxx[:, None]                 # ligands x receivers
    rss = np.maximum(syy[None, :] - b ** 2 * sxx[:, None], 1e-300)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = b / np.sqrt(rss / (n - 2) / sxx[:, None])
    t[sxx == 0] = np.nan
    return pd.DataFrame(t, index=W.columns, columns=profiles.columns)


def footprint_p(profiles: pd.DataFrame, top: int = TOP) -> pd.DataFrame:
    """ligands x receivers: the competitive footprint p-value (module docstring, use 2)."""
    a = activities(profiles, top)
    a = a[a.notna().all(1)]
    sd = a.std(1, ddof=0)
    z = a.sub(a.mean(1), axis=0).div(sd.where(sd > 0), axis=0).dropna()
    Z = z.to_numpy()
    n = Z.shape[0]
    P = np.empty_like(Z)
    for j in range(Z.shape[1]):
        col = np.sort(Z[:, j])
        ge = n - np.searchsorted(col, Z[:, j], side="left")      # #{L' : z(L') >= z(L)}, L included
        P[:, j] = ge / n
    return pd.DataFrame(P, index=z.index, columns=z.columns)
