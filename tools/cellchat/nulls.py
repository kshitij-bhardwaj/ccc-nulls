"""CellChat: Null A and Null B on CellChat's own score (communication probability ``prob``).

Null A is CellChat's own permutation test from ``computeCommunProb`` in CellChat 2.2.0.9001
(jinworks/CellChat@75253cd, R/modeling.R), vectorised. With ``rng="cellchat"`` (default) it draws
the permutations CellChat draws (R's ``set.seed`` + ``sample.int``, see ``cellchat_py.RRandom``),
so it returns the p-values CellChat reports. ``cellchat_py.compute_commun_prob`` is the literal,
loop-by-loop port it is checked against (``match_null_a.py``, ``test_adapter.py``):

    data      = data.signaling / max(data.signaling)            (all CellChatDB genes in the data)
    avg[g, k] = triMean of gene g in group k                   (collapse::fquantile, type 7)
    L[k]      = avg of the ligand; complexes: geometric mean of the subunits' avg
    R[k]      = same for the receptor, then * prod(1 + avg[co-activating receptors])
                                           / prod(1 + avg[co-inhibitory receptors])
    P1[s, t]  = L[s] R[t] / (Kh + L[s] R[t])                   (Hill function, Kh = 0.5, n = 1)
    P2[s, t]  = ag[s] ag[t],  ag = prod(1 + avg / (Kh + avg)) over agonist genes
    P3[s, t]  = an[s] an[t],  an = prod(Kh / (Kh + avg)) over antagonist genes
    prob      = P1 * P2 * P3 (* P4, population size, off by default)

    Null A: shuffle the group labels of all cells (B times), recompute prob, and
    p = #{b : prob_b > prob} / B, with p = 1 wherever prob = 0. No multiple-testing correction;
    CellChat calls an interaction significant when p < 0.05 (subsetCommunication).

Every floating-point step follows R's order of operations (R's two-pass mean, sequential
products, collapse's ``a + h * (b - a)`` interpolation) so that ties between permuted and
observed probabilities resolve the same way. tools/cellchat/README.md records the check.
"""
from __future__ import annotations


import numpy as np
import pandas as pd
from scipy import sparse

KH = 0.5
HILL_N = 1
SEP = ";"                                     # separator of gene lists in the adapter's extra columns


# ------------------------------------------------------------------ R arithmetic
def _r_mean(cols: list[np.ndarray]) -> np.ndarray:
    """R's mean() across ``cols``, elementwise: a sequential sum divided by n, then one
    correction pass, as in R's summary.c (long double == double on arm64 macOS)."""
    n = len(cols)
    s = cols[0].copy()
    for c in cols[1:]:
        s = s + c
    s = s / n
    with np.errstate(invalid="ignore"):
        t = cols[0] - s
        for c in cols[1:]:
            t = t + (c - s)
        return np.where(np.isfinite(s), s + t / n, s)


def _r_prod(cols: list[np.ndarray]) -> np.ndarray:
    """R's prod(): sequential product."""
    out = cols[0].copy()
    for c in cols[1:]:
        out = out * c
    return out


def _pow(x, n):
    return x if n == 1 else np.power(x, n)


# ------------------------------------------------------------------ trimean per group
class TriMean:
    """Exact per-group trimean (collapse::fquantile type 7, then R's mean of Q1, Q2, Q2, Q3) of
    the columns of a cells x genes sparse matrix, under any labelling of the cells.

    Values are non-negative, so a group's sorted values are (n - m) zeros followed by its m
    sorted non-zeros, and every order statistic is read off without sorting per group.
    """

    def __init__(self, U: sparse.spmatrix):
        U = sparse.csc_matrix(U, dtype=np.float64)
        U.eliminate_zeros()
        U.sort_indices()
        self.G = U.shape[1]
        col = np.repeat(np.arange(self.G), np.diff(U.indptr))
        order = np.lexsort((U.data, col))                  # by column, then by value
        self.col, self.row, self.val = col[order], U.indices[order], U.data[order]

    def __call__(self, codes: np.ndarray, K: int) -> np.ndarray:
        """K x G matrix of trimeans for the labelling ``codes`` (0..K-1 per cell)."""
        key = self.col * K + codes[self.row]
        vals = self.val[np.argsort(key, kind="stable")]    # value order kept within (gene, group)
        vals = np.append(vals, 0.0)                         # sentinel so empty genes index safely
        m = np.bincount(key, minlength=self.G * K).reshape(self.G, K)
        start = np.concatenate([[0], np.cumsum(m.ravel())[:-1]]).reshape(self.G, K)
        nk = np.bincount(codes, minlength=K).astype(np.int64)[None, :]
        zeros = nk - m

        def order_stat(k):                                  # k-th smallest value, 0-based
            nz = k >= zeros
            return np.where(nz, vals[np.where(nz, start + k - zeros, -1)], 0.0)

        q = []
        for Q in (0.25, 0.5, 0.75):
            h = (nk - 1) * Q                                 # collapse RETQSWITCH, type 7
            ih = h.astype(np.int64)
            frac = h - ih
            a = order_stat(np.broadcast_to(ih, m.shape))
            b = order_stat(np.broadcast_to(np.minimum(ih + 1, nk - 1), m.shape))
            last = np.broadcast_to((ih == nk - 1) | (frac <= 0.0), m.shape)
            q.append(np.where(last, a, a + frac * (b - a)))
        return _r_mean([q[0], q[1], q[1], q[2]]).T           # K x G


# ------------------------------------------------------------------ the CellChat score
def _genes(cell) -> list[str]:
    if cell is None or (isinstance(cell, float) and np.isnan(cell)) or cell == "":
        return []
    return str(cell).split(SEP)


GENE_COLS = ["ligand_genes", "receptor_genes", "agonist_genes", "antagonist_genes", "co_A_genes", "co_I_genes"]


def _lr_key(df: pd.DataFrame) -> pd.Series:
    if "interaction_name" in df.columns:
        return df["interaction_name"].astype(str)
    return df.ligand_complex.astype(str) + "|" + df.receptor_complex.astype(str)


def _lr_table(interactions: pd.DataFrame) -> pd.DataFrame:
    """One row per ligand-receptor pair with the gene lists CellChat's score needs.

    Rows from ``adapter.native_to_schema`` carry CellChat's own subunit and cofactor lists.
    Other rows (e.g. LIANA's) fall back to the complex names split on "_" and no cofactors.
    """
    df = interactions.assign(_key=_lr_key(interactions))
    if not all(c in df.columns for c in GENE_COLS):
        df = df.assign(ligand_genes=df.ligand_complex.astype(str).str.replace("_", SEP),
                       receptor_genes=df.receptor_complex.astype(str).str.replace("_", SEP),
                       agonist_genes="", antagonist_genes="", co_A_genes="", co_I_genes="")
    return df.drop_duplicates("_key").set_index("_key")[GENE_COLS].fillna("")


class CellChatScore:
    """CellChat's communication probability for many L-R pairs and all (source, target) groups."""

    def __init__(self, lr: pd.DataFrame, gene_pos: dict, Kh=KH, n=HILL_N, population_size=False):
        self.names = list(lr.index)
        self.Kh, self.n, self.population_size = Kh, n, population_size
        present = lambda gs: [gene_pos[g] for g in gs if g in gene_pos]
        self.L = [[gene_pos[g] for g in _genes(x)] for x in lr.ligand_genes]       # complexes: every
        self.R = [[gene_pos[g] for g in _genes(x)] for x in lr.receptor_genes]     # subunit
        self.coA = [present(_genes(x)) for x in lr.co_A_genes]                   # cofactors: those in
        self.coI = [present(_genes(x)) for x in lr.co_I_genes]                   # the data, as in R
        self.ag = [present(_genes(x)) for x in lr.agonist_genes]
        self.an = [present(_genes(x)) for x in lr.antagonist_genes]

    @staticmethod
    def _unit(avg, idx):
        """A ligand's / receptor's expression: its gene's avg, or the subunits' geometric mean."""
        if len(idx) == 1:
            return avg[:, idx[0]]
        with np.errstate(divide="ignore"):
            return np.exp(_r_mean([np.log(avg[:, g]) for g in idx]))

    def __call__(self, avg: np.ndarray, sizes: np.ndarray | None = None) -> np.ndarray:
        """nLR x K x K probabilities from the K x G group averages ``avg``."""
        K, Kh, n = avg.shape[0], self.Kh, self.n
        Khn = Kh ** n
        out = np.empty((len(self.names), K, K))
        for i in range(len(self.names)):
            Lv = self._unit(avg, self.L[i])
            Rv = self._unit(avg, self.R[i])
            coA = _r_prod([1 + avg[:, g] for g in self.coA[i]]) if self.coA[i] else 1.0
            coI = _r_prod([1 + avg[:, g] for g in self.coI[i]]) if self.coI[i] else 1.0
            Rv = Rv * coA / coI
            LR = Lv[:, None] * Rv[None, :]
            P = _pow(LR, n) / (Khn + _pow(LR, n))
            if self.ag[i]:
                a = _r_prod([1 + _pow(avg[:, g], n) / (Khn + _pow(avg[:, g], n)) for g in self.ag[i]])
                P = P * (a[:, None] * a[None, :])
            if self.an[i]:
                a = _r_prod([Khn / (Khn + _pow(avg[:, g], n)) for g in self.an[i]])
                P = P * (a[:, None] * a[None, :])
            if self.population_size:
                f = sizes / sizes.sum()
                P = P * (f[:, None] * f[None, :])
            out[i] = P
        return out


# ------------------------------------------------------------------ Null A
def signaling_genes() -> list[str]:
    """Genes of ``data.signaling``: every gene of the default CellChatDB (``extractGene``)."""
    from tools.cellchat import cellchat_py as cc

    return cc.extract_gene(cc.subset_db(cc.load_db("human")))


def null_a(adata, interactions: pd.DataFrame, groupby: str = "label", *, B: int = 1000, seed: int = 0,
           rng: str = "cellchat", permutations: np.ndarray | None = None,
           genes_signaling: list[str] | None = None, population_size: bool = False,
           Kh: float = KH, n: int = HILL_N) -> pd.DataFrame:
    """Null A: CellChat's permutation test on CellChat's score.

    Parameters
    ----------
    adata
        ``X`` = the log-normalised expression CellChat scored; ``obs[groupby]`` = cell groups.
    interactions
        Common-schema rows, ideally from ``adapter.native_to_schema`` (they carry CellChat's
        subunit and cofactor lists). Rows of the same L-R pair share one computation.
    B, seed
        Number of permutations and their seed. CellChat's defaults are B = 100, seed = 1.
    rng
        ``"cellchat"``: CellChat's permutations, ``set.seed(seed); replicate(B, sample.int(nC, nC))``
        with R's generator. ``"numpy"``: independent permutations from numpy.
    permutations
        Optional cells x B matrix of 1-based cell indices; overrides ``B``, ``seed`` and ``rng``.
    genes_signaling
        Genes of ``data.signaling`` (all CellChatDB genes in the data); their maximum scales the
        data. Default: the genes of the default CellChatDB.

    Returns
    -------
    ``interactions`` plus
        ``score_null_a``       CellChat's prob recomputed here
        ``n_exceed_null_a``    #{b : prob_b > prob}
        ``p_null_a_cellchat``  CellChat's p-value: n_exceed / B, and 1 where prob = 0
        ``p_null_a``           the repository's convention, (1 + n_exceed) / (B + 1), 1 where prob = 0
    """
    labels = adata.obs[groupby].astype(str).to_numpy()
    levels = np.unique(labels)
    codes = np.searchsorted(levels, labels)
    K = len(levels)
    genes = pd.Index(adata.var_names.astype(str))

    lr = _lr_table(interactions)
    sig = genes_signaling if genes_signaling is not None else signaling_genes()
    sig_idx = genes.get_indexer(pd.Index(sig).unique())
    X = sparse.csc_matrix(adata.X)
    data_max = float(X[:, sig_idx[sig_idx >= 0]].max())

    need = sorted({g for c in GENE_COLS for x in lr[c] for g in _genes(x)} & set(genes))
    missing = {g for c in ("ligand_genes", "receptor_genes") for x in lr[c] for g in _genes(x)} - set(need)
    if missing:
        raise ValueError(f"ligand / receptor genes not in the data: {sorted(missing)[:10]}")
    U = sparse.csc_matrix(X[:, genes.get_indexer(need)], dtype=np.float64)
    U.data = U.data / data_max          # true division, as R's data / max(data); scipy's U / x multiplies by 1/x
    tm = TriMean(U)
    score = CellChatScore(lr, {g: j for j, g in enumerate(need)}, Kh=Kh, n=n, population_size=population_size)
    sizes = np.bincount(codes, minlength=K).astype(float)

    P_obs = score(tm(codes, K), sizes)
    exceed = np.zeros(P_obs.shape, dtype=np.int64)
    if permutations is None and rng == "cellchat":
        from tools.cellchat.cellchat_py import r_permutations

        permutations = r_permutations(len(codes), B, seed)
    if permutations is not None:
        permutations = np.asarray(permutations)
        B = permutations.shape[1]
    gen = np.random.default_rng(seed)
    for b in range(B):
        cb = codes[permutations[:, b] - 1] if permutations is not None else gen.permutation(codes)
        exceed += score(tm(cb, K), sizes) > P_obs

    out = interactions.copy()
    i = pd.Index(score.names).get_indexer(_lr_key(out))
    s = np.searchsorted(levels, out.source.astype(str).to_numpy())
    t = np.searchsorted(levels, out.target.astype(str).to_numpy())
    prob, k = P_obs[i, s, t], exceed[i, s, t]
    zero = prob == 0
    out["score_null_a"] = prob
    out["n_exceed_null_a"] = k
    out["p_null_a_cellchat"] = np.where(zero, 1.0, k / B)
    out["p_null_a"] = np.where(zero, 1.0, (1 + k) / (B + 1))
    out.attrs["null_a"] = dict(B=int(B), seed=int(seed), rng=rng if permutations is None or rng == "cellchat" else "given",
                               data_max=data_max, n_lr=len(score.names), groups=[str(x) for x in levels])
    return out


def null_b(adata, interactions: pd.DataFrame, groupby: str = "label", *, B: int = 1000,
           seed: int = 0) -> pd.DataFrame:
    """Null B: keep the labels and rank the observed score against abundance-matched decoy genes.

    Returns ``interactions`` with ``p_null_b`` added.
    """
    raise NotImplementedError("cellchat: null_b is not implemented yet (see tools/cellchat/nulls.py)")
