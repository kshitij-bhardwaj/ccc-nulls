"""Null B as a composite null: is the sender -> receiver edge specific on both sides? Any tool.

A call "s signals to t through L-R" needs two things: the ligand side is specific to s AND the
receptor side is specific to t. Null A shuffles the labels of both sides at once, so it rejects as
soon as either side is specific (a ligand made by s and a receptor made everywhere is called for
every receiver). Null B tests the null of the two-part claim, a union:

    H0 = H0_L or H0_R
    H0_L: the sender side's group profile is unrelated to the labels (receiver side as observed)
    H0_R: the receiver side's group profile is unrelated to the labels (sender side as observed)

p_L shuffles the labels for the sender side only and recomputes the tool's score, p_R the same for
the receiver side, and p = max(p_L, p_R): the intersection-union test, level alpha for the union
without adjustment (Berger 1982, Technometrics 24:295; Berger & Hsu 1996, Statist. Sci. 11:283). It
is the partial-conjunction p-value for "2 of 2 components non-null" (Benjamini & Heller 2008,
Biometrics 64:1215), and the joint-significance test of no mediation (MacKinnon et al. 2002,
Psychol. Methods 7:83). No gene bins, no decoys: the tool's own score is tested.

How a tool plugs in: its adapter builds a :class:`SideScore` from two callables::

    summarize(codes, n_types) -> per-group statistics of the genes, for the labelling ``codes``
    score(sender, receiver)   -> the tool's score, with the sender side (ligand, sender factors)
                                 read from ``sender`` and the receiver side from ``receiver``

such that ``score(summarize(c), summarize(c))`` is the tool's score. The engine never looks inside
the statistics, so CellChat passes trimeans and its Hill function, CellPhoneDB means and their
average. Scores of the form ``f(M, li, ri, si, ti)`` on cluster means (:mod:`ccc_nulls.scores`)
get one for free from :func:`mean_side_score`, and :func:`null_b_table` runs the whole thing on
common-schema candidates.

Null B methods are registered by name (:func:`register`). A method takes ``(side, codes, n_types,
**options)``, where adapters pass ``B, seed, permutations, delta, ties`` (a method ignores what does
not apply to it), and returns a dict with at least ``score`` and ``p`` (see :func:`to_columns`). Adapters call
``get(method)(side, ...)``, so a Null B added here runs on every tool that provides a side score,
with no change to the tool's folder.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
import pandas as pd

from ._matrix import as_matrix, cluster_means, columns
from .scores import get_score
from .stats import bh

SUFFIX = "null_b"


@dataclass(frozen=True)
class SideScore:
    """A tool's score split into its sender and receiver sides (see module docstring)."""

    summarize: Callable[[np.ndarray, int], Any]
    score: Callable[[Any, Any], np.ndarray]
    name: str = ""


# ------------------------------------------------------------------ registry of Null B methods
METHODS: dict[str, Callable] = {}


def register(name: str):
    """Decorator: make a Null B method available to every adapter as ``get(name)``."""
    def deco(fn):
        if name in METHODS:
            raise ValueError(f"Null B method {name!r} is already registered")
        METHODS[name] = fn
        return fn
    return deco


def get(method: str | Callable = "composite") -> Callable:
    """The Null B method registered as ``method``, or ``method`` itself if it is a callable."""
    if callable(method):
        return method
    try:
        return METHODS[method]
    except KeyError:
        raise ValueError(f"unknown Null B method {method!r}; registered: {sorted(METHODS)}") from None


# ------------------------------------------------------------------ the composite null
@register("composite")
def composite(side: SideScore, codes, n_types: int, *, B: int = 1000, seed: int = 0,
              permutations: np.ndarray | None = None, delta: float | tuple = 0.0,
              ties: str = "count") -> dict:
    """p_L, p_R and p = max(p_L, p_R) for every entry of ``side.score``'s output.

    Parameters
    ----------
    side
        The tool's :class:`SideScore`. Its score is non-negative; 0 means no signal, and p = 1 there.
    codes, n_types
        Group code (0..n_types-1) of every cell.
    B, seed
        Label permutations, drawn with ``numpy.random.default_rng(seed)``. Both sides use the same
        permutations, as does Null A for equal ``B`` and ``seed``.
    permutations
        Optional cells x B matrix of 0-based cell indices; overrides ``B`` and ``seed``.
    delta
        Interval null: each side may add at most delta to the score over random labels
        (minimum-effect testing; Hodges & Lehmann 1954; McCarthy & Smyth 2009). A tuple runs several
        deltas on one set of permutations. Sensitivity setting; 0 is the null above.
    ties
        ``"count"``: p = (1 + #{b : s_b >= s - delta}) / (B + 1), the repository's convention.
        ``"ignore"``: count s_b > s - delta only, as CellChat's own test does.

    Returns
    -------
    dict of arrays shaped like ``side.score``'s output: ``score``, ``p_L``, ``p_R``, ``p``,
    ``limiting_side`` ("L" or "R": the side with the larger p), ``n_pos_L`` / ``n_pos_R``
    (#{b : s_b > 0} with that side permuted, for :func:`testable`), and ``p_delta`` (delta -> p);
    plus ``B`` and ``deltas``.
    """
    if ties not in ("count", "ignore"):
        raise ValueError("ties must be 'count' or 'ignore'")
    beats = np.greater_equal if ties == "count" else np.greater
    deltas = (delta,) if np.isscalar(delta) else tuple(delta)
    codes = np.asarray(codes)
    obs = side.summarize(codes, n_types)
    s_obs = np.asarray(side.score(obs, obs))
    kL = {d: np.zeros(s_obs.shape, dtype=np.int64) for d in deltas}
    kR = {d: np.zeros(s_obs.shape, dtype=np.int64) for d in deltas}
    nL, nR = np.zeros(s_obs.shape, dtype=np.int64), np.zeros(s_obs.shape, dtype=np.int64)
    if permutations is not None:
        permutations = np.asarray(permutations)
        B = permutations.shape[1]
    gen = np.random.default_rng(seed)
    for b in range(B):
        perm = side.summarize(codes[permutations[:, b]] if permutations is not None else gen.permutation(codes), n_types)
        s_L = np.asarray(side.score(perm, obs))              # sender side shuffled
        s_R = np.asarray(side.score(obs, perm))              # receiver side shuffled
        nL += s_L > 0
        nR += s_R > 0
        for d in deltas:
            kL[d] += beats(s_L, s_obs - d)
            kR[d] += beats(s_R, s_obs - d)

    pval = lambda k: np.where(s_obs == 0, 1.0, (1 + k) / (B + 1))
    p_delta = {d: np.maximum(pval(kL[d]), pval(kR[d])) for d in deltas}
    pL, pR = pval(kL[deltas[0]]), pval(kR[deltas[0]])
    return dict(score=s_obs, p_L=pL, p_R=pR, p=p_delta[deltas[0]], limiting_side=np.where(pL >= pR, "L", "R"),
                n_pos_L=nL, n_pos_R=nR, p_delta=p_delta, B=int(B), deltas=[float(d) for d in deltas])


def to_columns(res: dict, index=slice(None)) -> dict:
    """The result as output columns (``score_null_b``, ``p_null_b``, ...), read at ``index``
    (e.g. the (pair, source, target) indices of each output row into the score's array)."""
    out = {f"score_{SUFFIX}": res["score"][index]}
    if "n_pos_L" in res:
        out[f"n_pos_L_{SUFFIX}"], out[f"n_pos_R_{SUFFIX}"] = res["n_pos_L"][index], res["n_pos_R"][index]
    if "p_L" in res:
        out[f"p_L_{SUFFIX}"], out[f"p_R_{SUFFIX}"] = res["p_L"][index], res["p_R"][index]
    out[f"p_{SUFFIX}"] = res["p"][index]
    if "limiting_side" in res:
        out[f"limiting_side_{SUFFIX}"] = res["limiting_side"][index]
    for d, p in list(res.get("p_delta", {}).items())[1:]:
        out[f"p_{SUFFIX}_delta={d:g}"] = p[index]
    return out


def testable(score, n_pos_L, n_pos_R, *more) -> np.ndarray:
    """The rows to correct over: some labelling (observed, or either side permuted; ``more``: other
    nulls' permutations, e.g. Null A's #{s_b > 0} > 0) gives score > 0.

    Not the rows the tool reports (score > 0 under the observed labels): that set needs the receiver
    side to be non-zero in the observed labelling, so it selects on the side H0_R calls random, and
    p_R is then too small. Membership of the testable set does not depend on which labelling is the
    observed one, and its rows with score = 0 have p = 1, which only makes BH more conservative.
    """
    out = (np.asarray(score) > 0) | (np.asarray(n_pos_L) > 0) | (np.asarray(n_pos_R) > 0)
    for m in more:
        out = out | np.asarray(m, bool)
    return out


# ------------------------------------------------------------------ scores on cluster means
def mean_side_score(X, li, ri, si, ti, score="mean") -> SideScore:
    """Side score for ``f(M, li, ri, si, ti)`` (:mod:`ccc_nulls.scores`, or a callable that reads
    only ligand columns at source rows and receptor columns at target rows). Output: one value per
    candidate. The sender and receiver cluster means are stacked side by side, so ``f`` is used as is.
    """
    f = get_score(score)
    used, inv = np.unique(np.concatenate([li, ri]), return_inverse=True)
    li_, ri_ = inv[: len(li)], inv[len(li):] + len(used)
    Xu = columns(as_matrix(X), used)
    return SideScore(summarize=lambda codes, K: cluster_means(Xu, codes, K),
                     score=lambda Ms, Mr: f(np.hstack([Ms, Mr]), li_, ri_, si, ti),
                     name=getattr(f, "__name__", "score"))


def null_b_table(X, genes, labels, candidates: pd.DataFrame, *, score="mean", method: str | Callable = "composite",
                 alpha: float = 0.05, **options) -> pd.DataFrame:
    """Null B on common-schema candidates (``ligand, receptor, source, target``) for a score on
    cluster means, the entry point for tools without their own side score.

    Returns the candidates that map onto X plus the method's columns (:func:`to_columns`),
    ``q_null_b`` (BH over :func:`testable` rows when the method reports them, else over all rows; 1
    outside) and ``call_null_b`` (q <= alpha). ``options`` go to the method (B, seed, delta, ...).
    """
    X = as_matrix(X)
    genes = pd.Index(np.asarray(genes).astype(str))
    labels = np.asarray(labels).astype(str)
    if X.shape != (len(labels), len(genes)):
        raise ValueError(f"X is {X.shape} but there are {len(labels)} labels and {len(genes)} genes")
    types = pd.Index(sorted(set(labels)))
    c = candidates.copy()
    idx = [genes.get_indexer(c.ligand.astype(str)), genes.get_indexer(c.receptor.astype(str)),
           types.get_indexer(c.source.astype(str)), types.get_indexer(c.target.astype(str))]
    ok = np.logical_and.reduce([i >= 0 for i in idx])
    if not ok.all():
        warnings.warn(f"dropping {int((~ok).sum())} candidates whose genes or cell types are not in X", stacklevel=2)
    c, (li, ri, si, ti) = c.loc[ok].copy(), [i[ok] for i in idx]
    if len(c) == 0:
        raise ValueError("no candidate maps onto X")

    res = get(method)(mean_side_score(X, li, ri, si, ti, score), types.get_indexer(labels), len(types), **options)
    for k, v in to_columns(res).items():
        c[k] = v
    fam = testable(res["score"], res["n_pos_L"], res["n_pos_R"]) if "n_pos_L" in res else np.ones(len(c), bool)
    q = np.ones(len(c))
    q[fam] = bh(res["p"][fam])
    c[f"q_{SUFFIX}"], c[f"call_{SUFFIX}"] = q, q <= alpha
    c.attrs[SUFFIX] = dict(method=method if isinstance(method, str) else getattr(method, "__name__", "callable"),
                           B=res.get("B"), deltas=res.get("deltas"), n_testable=int(fam.sum()))
    return c
