"""Simulator with known ground truth, for calibrating null models.

Regimes
-------
S0  ``structured=False``, labels shuffled afterwards: no cell-type structure, no CCC.
S1  ``structured=True, planted=False``: every gene carries a cell-type effect
    ``delta ~ N(0, tau^2)`` (log scale), but no LR pair is special. Every call is a
    false positive.
S2  ``structured=True, planted=True``: S1 plus ``n_planted`` events with the ligand up
    in X and the receptor up in Y by ``plant_lfc``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class SimData:
    X: np.ndarray                     # cells x genes, log1p(CP10k)
    genes: np.ndarray
    labels: np.ndarray                # cell-type names
    pairs: pd.DataFrame               # the LR "resource": ligand, receptor
    truth: set = field(default_factory=set)      # {(ligand, receptor, source, target)}
    lig_up: set = field(default_factory=set)     # {(ligand, source)} planted
    rec_up: set = field(default_factory=set)     # {(receptor, target)} planted


def simulate(seed=0, n_genes=8000, sizes=(500, 400, 300, 200, 150, 100), tau=0.25,
             n_pairs=400, n_planted=60, plant_lfc=1.0, dispersion=0.3,
             structured=True, planted=True, shuffle_labels=False) -> SimData:
    """Negative-binomial counts with a broad abundance spectrum, CP10k + log1p normalised."""
    rng = np.random.default_rng(seed)
    K, n = len(sizes), sum(sizes)
    codes = np.repeat(np.arange(K), sizes)
    log_mu = rng.normal(-1.0, 1.5, n_genes)
    delta = rng.normal(0, tau, (K, n_genes)) if structured else np.zeros((K, n_genes))
    pairs = np.array([rng.choice(n_genes, 2, replace=False) for _ in range(n_pairs)])
    truth, lig_up, rec_up = set(), set(), set()
    genes = np.array([f"G{i:05d}" for i in range(n_genes)])
    types = np.array([f"T{k}" for k in range(K)])
    if planted:
        for j in rng.choice(n_pairs, n_planted, replace=False):
            x, y = rng.choice(K, 2, replace=True)
            L, R = pairs[j]
            delta[x, L] += plant_lfc
            delta[y, R] += plant_lfc
            truth.add((genes[L], genes[R], types[x], types[y]))
            lig_up.add((genes[L], types[x]))
            rec_up.add((genes[R], types[y]))
    size_f = rng.lognormal(0, 0.3, n)
    mu = np.exp(log_mu[None, :] + delta[codes]) * size_f[:, None]
    counts = rng.poisson(rng.gamma(1 / dispersion, dispersion * mu))
    lib = np.maximum(counts.sum(1, keepdims=True), 1)
    X = np.log1p(counts / lib * 1e4)
    if shuffle_labels:
        codes = np.random.default_rng(seed + 99).permutation(codes)
    return SimData(X, genes, types[codes], pd.DataFrame({"ligand": genes[pairs[:, 0]], "receptor": genes[pairs[:, 1]]}),
                   truth, lig_up, rec_up)


def candidates_from_pairs(X, genes, labels, pairs: pd.DataFrame, min_frac=0.10) -> pd.DataFrame:
    """Every (L, R, X, Y) passing a CellPhoneDB-style expression filter, in schema form."""
    types = np.array(sorted(set(labels)))
    gi = pd.Index(genes)
    frac = np.stack([(X[labels == t] > 0).mean(0) for t in types])
    rows = []
    for L, R in zip(gi.get_indexer(pairs.ligand), gi.get_indexer(pairs.receptor)):
        for x in np.where(frac[:, L] >= min_frac)[0]:
            for y in np.where(frac[:, R] >= min_frac)[0]:
                rows.append((genes[L], genes[R], types[x], types[y]))
    return pd.DataFrame(rows, columns=["ligand", "receptor", "source", "target"])


def label_truth(cands: pd.DataFrame, sim: SimData) -> np.ndarray:
    """'true', 'partial' (one side planted), or 'null' for each candidate."""
    key = list(zip(cands.ligand, cands.receptor, cands.source, cands.target))
    true = np.array([k in sim.truth for k in key])
    part = np.array([(l, s) in sim.lig_up or (r, t) in sim.rec_up for l, r, s, t in key]) & ~true
    return np.where(true, "true", np.where(part, "partial", "null"))
