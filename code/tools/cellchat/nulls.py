"""CellChat: Null A and Null B on CellChat's own score (communication probability ``prob``).

Null A is CellChat's own permutation test from ``computeCommunProb`` in CellChat 2.1.2
(jinworks/CellChat, release v2.1.2, R/modeling.R), vectorised. With ``rng="cellchat"`` (default) it draws
the permutations CellChat draws (R's ``set.seed`` + ``sample.int``, see ``cellchat_py.RRandom``),
so it returns the p-values CellChat reports. ``cellchat_py.compute_commun_prob`` is the literal,
loop-by-loop port it is checked against (``match_null_a.py``, ``test_adapter.py``):

    data      = data.signaling / max(data.signaling)            (all CellChatDB genes in the data)
    avg[g, k] = triMean of gene g in group k                   (stats::quantile, type 7)
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

    Null B (``adapter.null_b``): the composite null "the sender side is not specific to s OR the
    receiver side is not specific to t", run by the shared engine ``ccc_nulls.composite`` on
    CellChat's score split into sides (``CellChatScore(avg, avg_receiver=...)``, wrapped by
    ``adapter.side_score``). p = max(p_L, p_R) (intersection-union test). BH over
    ``adapter.testable`` rows.

Every floating-point step follows R's order of operations (R's two-pass mean, sequential
products, R's quantile interpolation) so that ties between permuted and
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
    """Exact per-group trimean (type-7 quartiles, then R's mean of Q1, Q2, Q2, Q3) of the columns
    of a cells x genes sparse matrix, under any labelling of the cells. ``quantile="collapse"``
    interpolates as collapse::fquantile (CellChat >= 2), ``"stats"`` as stats::quantile (1.x).

    Values are non-negative, so a group's sorted values are (n - m) zeros followed by its m
    sorted non-zeros, and every order statistic is read off without sorting per group.
    """

    def __init__(self, U: sparse.spmatrix, quantile: str = "collapse"):
        self.quantile = quantile
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
            if self.quantile == "collapse":
                q.append(np.where(last, a, a + frac * (b - a)))
            else:                                            # stats::quantile, only where b != a
                q.append(np.where(last | (b == a), a, (1 - frac) * a + frac * b))
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

    def __call__(self, avg: np.ndarray, sizes: np.ndarray | None = None, avg_receiver: np.ndarray | None = None) -> np.ndarray:
        """nLR x K x K probabilities from the K x G group averages ``avg``. With ``avg_receiver``,
        the receiver side (receptor, co-receptors, the receiver's agonist / antagonist factor) is
        read from it and the sender side (ligand, the sender's factor) from ``avg`` (Null B)."""
        a_s, a_r = avg, avg if avg_receiver is None else avg_receiver
        K, Kh, n = avg.shape[0], self.Kh, self.n
        Khn = Kh ** n
        out = np.empty((len(self.names), K, K))
        for i in range(len(self.names)):
            Lv = self._unit(a_s, self.L[i])
            Rv = self._unit(a_r, self.R[i])
            coA = _r_prod([1 + a_r[:, g] for g in self.coA[i]]) if self.coA[i] else 1.0
            coI = _r_prod([1 + a_r[:, g] for g in self.coI[i]]) if self.coI[i] else 1.0
            Rv = Rv * coA / coI
            LR = Lv[:, None] * Rv[None, :]
            P = _pow(LR, n) / (Khn + _pow(LR, n))
            if self.ag[i]:
                ag = lambda a: _r_prod([1 + _pow(a[:, g], n) / (Khn + _pow(a[:, g], n)) for g in self.ag[i]])
                P = P * (ag(a_s)[:, None] * ag(a_r)[None, :])
            if self.an[i]:
                an = lambda a: _r_prod([Khn / (Khn + _pow(a[:, g], n)) for g in self.an[i]])
                P = P * (an(a_s)[:, None] * an(a_r)[None, :])
            if self.population_size:
                f = sizes / sizes.sum()
                P = P * (f[:, None] * f[None, :])
            out[i] = P
        return out


class CellChatScoreV1:
    """CellChat 1.0.0's probability: ligand and receptor expression are formed per cell (complex =
    geometric mean of its subunits in each cell; receptor * prod(1 + co-A) / prod(1 + co-I) in each
    cell), then averaged per group; agonists and antagonists use per-gene group averages."""

    def __init__(self, lr: pd.DataFrame, U: sparse.csc_matrix, gene_pos: dict, Kh=KH):
        self.names, self.Kh = list(lr.index), Kh
        present = lambda gs: [gene_pos[g] for g in gs if g in gene_pos]
        col = lambda g: U[:, g].toarray().ravel()
        units, self.Li, self.Ri = {}, [], []
        for lig, rec, coA, coI in zip(lr.ligand_genes, lr.receptor_genes, lr.co_A_genes, lr.co_I_genes):
            for key, side in [(("L", lig), self.Li), (("R", rec, coA, coI), self.Ri)]:
                if key not in units:
                    units[key] = len(units)
                side.append(units[key])
        vecs = []
        for key in units:
            idx = [gene_pos[g] for g in _genes(key[1])]
            with np.errstate(divide="ignore"):
                v = col(idx[0]) if len(idx) == 1 else np.exp(_r_mean([np.log(col(g)) for g in idx]))
            if key[0] == "R":
                a, i = present(_genes(key[2])), present(_genes(key[3]))
                v = v * (_r_prod([1 + col(g) for g in a]) if a else 1.0) / (_r_prod([1 + col(g) for g in i]) if i else 1.0)
            vecs.append(sparse.csc_matrix(v[:, None]))
        self.units = sparse.hstack(vecs, format="csc") if vecs else sparse.csc_matrix((U.shape[0], 0))
        self.ag = [present(_genes(x)) for x in lr.agonist_genes]
        self.an = [present(_genes(x)) for x in lr.antagonist_genes]

    def __call__(self, avg_units: np.ndarray, avg_genes: np.ndarray) -> np.ndarray:
        K, Kh = avg_units.shape[0], self.Kh
        out = np.empty((len(self.names), K, K))
        for i in range(len(self.names)):
            LR = avg_units[:, self.Li[i]][:, None] * avg_units[:, self.Ri[i]][None, :]
            P = LR / (Kh + LR)
            if self.ag[i]:
                a = _r_prod([1 + avg_genes[:, g] / (Kh + avg_genes[:, g]) for g in self.ag[i]])
                P = P * (a[:, None] * a[None, :])
            if self.an[i]:
                a = _r_prod([Kh / (Kh + avg_genes[:, g]) for g in self.an[i]])
                P = P * (a[:, None] * a[None, :])
            out[i] = P
        return out


# ------------------------------------------------------------------ Null A
def signaling_genes() -> list[str]:
    """Genes of ``data.signaling``: every gene of the default CellChatDB (``extractGene``)."""
    from tools.cellchat import cellchat_py as cc

    return cc.extract_gene(cc.subset_db(cc.load_db("human")))


def _setup(adata, interactions, groupby, genes_signaling):
    """What both nulls start from: group codes, the L-R table, and the data CellChat scores
    (the needed genes of data.signaling divided by its maximum)."""
    labels = adata.obs[groupby].astype(str).to_numpy()
    levels = np.unique(labels)
    codes = np.searchsorted(levels, labels)
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
    gene_pos = {g: j for j, g in enumerate(need)}
    return levels, codes, len(levels), lr, U, gene_pos, data_max


def _rows(out, names, levels):
    """Indices (L-R pair, source, target) of each row of ``out`` into nLR x K x K arrays."""
    i = pd.Index(names).get_indexer(_lr_key(out))
    s = np.searchsorted(levels, out.source.astype(str).to_numpy())
    t = np.searchsorted(levels, out.target.astype(str).to_numpy())
    return i, s, t


def null_a(adata, interactions: pd.DataFrame, groupby: str = "label", *, B: int = 1000, seed: int = 0,
           rng: str = "cellchat", permutations: np.ndarray | None = None,
           genes_signaling: list[str] | None = None, population_size: bool = False,
           Kh: float = KH, n: int = HILL_N, version: str = "2.1.2", keep_top: int = 0) -> pd.DataFrame:
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
    version
        ``"2.1.2"`` (default): the latest CellChat release; p counts permutations with prob_b > prob.
        ``"2.2.0.9001"``: GitHub main; as 2.1.2 but triMean uses collapse::fquantile.
        ``"1.0.0"``: CellChat 1.0.0, the 2021 paper's version: complexes and co-receptors per cell,
        stats::quantile, and p counts prob_b >= prob (ties count against the interaction).
    keep_top
        If > 0, also keep the ``keep_top`` largest permuted probabilities of every interaction with
        prob > 0 (``perm_top_null_a``). Outputs are unchanged.

    Returns
    -------
    ``interactions`` plus ``perm_mean_null_a``, ``perm_sd_null_a`` and ``z_null_a``
    (``(prob - perm_mean) / perm_sd``, NaN where every permutation gives the same prob), the
    permutation counts and moments Null B's p-values use (``n_tie_null_a`` = #{prob_b = prob},
    ``perm_npos_null_a`` = #{prob_b > 0}, ``perm_pos_{mean,sd,skew}_null_a`` of the positive prob_b,
    ``perm_top_null_a`` = the ``keep_top`` largest prob_b, descending), and
        ``score_null_a``       CellChat's prob recomputed here
        ``n_exceed_null_a``    #{b : prob_b > prob}
        ``p_null_a_cellchat``  CellChat's p-value: n_exceed / B, and 1 where prob = 0
        ``p_null_a``           the repository's convention, (1 + n_exceed) / (B + 1), 1 where prob = 0
    """
    levels, codes, K, lr, U, gene_pos, data_max = _setup(adata, interactions, groupby, genes_signaling)
    sizes = np.bincount(codes, minlength=K).astype(float)
    if version in ("2.1.2", "2.2.0.9001"):
        tm = TriMean(U, quantile="stats" if version == "2.1.2" else "collapse")
        sc = CellChatScore(lr, gene_pos, Kh=Kh, n=n, population_size=population_size)
        score_of = lambda cd: sc(tm(cd, K), sizes)
        beats = np.greater
    elif version == "1.0.0":
        if population_size or n != 1:
            raise NotImplementedError("version 1.0.0 is implemented for population_size=False and n=1")
        sc = CellChatScoreV1(lr, U, gene_pos, Kh=Kh)
        tm_units, tm_genes = TriMean(sc.units, quantile="stats"), TriMean(U, quantile="stats")
        score_of = lambda cd: sc(tm_units(cd, K), tm_genes(cd, K))
        beats = np.greater_equal
    else:
        raise ValueError("version must be '2.1.2', '2.2.0.9001' or '1.0.0'")
    score = sc

    P_obs = score_of(codes)
    exceed = np.zeros(P_obs.shape, dtype=np.int64)
    tie, npos = np.zeros(P_obs.shape, dtype=np.int64), np.zeros(P_obs.shape, dtype=np.int64)
    s1, s2, s3 = np.zeros(P_obs.shape), np.zeros(P_obs.shape), np.zeros(P_obs.shape)  # sums of prob_b^1,2,3
    scored = np.flatnonzero(P_obs > 0) if keep_top else np.array([], dtype=np.int64)
    top, buf = np.empty((len(scored), 0)), []                 # running top-k of prob_b, scored cells only
    if permutations is None and rng == "cellchat":
        from tools.cellchat.cellchat_py import r_permutations

        permutations = r_permutations(len(codes), B, seed)
    if permutations is not None:
        permutations = np.asarray(permutations)
        B = permutations.shape[1]
    if not keep_top:                                       # permutations on $CCC_NULLS_WORKERS processes
        from ccc_nulls import composite

        def counts(lo, hi):
            ex, ti, npo = (np.zeros(P_obs.shape, dtype=np.int64) for _ in range(3))
            a1, a2, a3 = np.zeros(P_obs.shape), np.zeros(P_obs.shape), np.zeros(P_obs.shape)
            g = np.random.default_rng(seed)
            for b in range(hi):
                if permutations is not None:
                    if b < lo:
                        continue
                    cb = codes[permutations[:, b] - 1]
                else:
                    cb = g.permutation(codes)
                    if b < lo:
                        continue
                Pb = score_of(cb)
                ex += beats(Pb - P_obs, 0)
                ti += Pb == P_obs
                npo += Pb > 0
                a1 += Pb
                a2 += Pb * Pb
                a3 += Pb ** 3
            return ex, ti, npo, a1, a2, a3

        exceed, tie, npos, s1, s2, s3 = composite.map_permutations(counts, int(B), composite.n_workers())
        B_loop = 0
    else:
        B_loop = B
    gen = np.random.default_rng(seed)
    for b in range(B_loop):
        cb = codes[permutations[:, b] - 1] if permutations is not None else gen.permutation(codes)
        Pb = score_of(cb)
        exceed += beats(Pb - P_obs, 0)                    # R: rowSums(Pboot - Pnull > 0), or >= 0 in 1.0
        tie += Pb == P_obs
        npos += Pb > 0
        s1 += Pb
        s2 += Pb * Pb
        s3 += Pb ** 3                                     # zeros add nothing: these are also the sums over prob_b > 0
        if keep_top:
            buf.append(Pb.ravel()[scored])
            if len(buf) == keep_top or b == B - 1:        # merge in blocks: partition keeps the k largest
                top = np.hstack([top, np.stack(buf, 1)])
                top = -np.partition(-top, min(keep_top, top.shape[1]) - 1, axis=1)[:, :keep_top]
                buf = []

    out = interactions.copy()
    i, s, t = _rows(out, score.names, levels)
    prob, k = P_obs[i, s, t], exceed[i, s, t]
    zero = prob == 0
    out["score_null_a"] = prob
    out["n_exceed_null_a"] = k
    out["p_null_a_cellchat"] = np.where(zero, 1.0, k / B)
    out["p_null_a"] = np.where(zero, 1.0, (1 + k) / (B + 1))
    # permutation moments of prob, for the standardized statistic z_null_a
    mean = s1 / B
    sd = np.sqrt(np.maximum(s2 - B * mean * mean, 0) / max(B - 1, 1))
    m_, sd_ = mean[i, s, t], sd[i, s, t]
    out["perm_mean_null_a"] = m_
    out["perm_sd_null_a"] = sd_
    with np.errstate(divide="ignore", invalid="ignore"):
        out["z_null_a"] = np.where(sd_ > 0, (prob - m_) / sd_, np.nan)    # NaN where the permutation null is degenerate
        # moments of the positive part of the permutation distribution (Null B, Pearson type III tail)
        c = np.maximum(npos, 1)
        pm, p2, p3 = s1 / c, s2 / c, s3 / c
        pv = np.maximum(p2 - pm * pm, 0)
        psd = np.sqrt(pv * c / np.maximum(c - 1, 1))
        pskew = np.where(pv > 0, (p3 - 3 * pm * p2 + 2 * pm ** 3) / pv ** 1.5, 0.0)
    out["n_tie_null_a"] = tie[i, s, t]
    out["perm_npos_null_a"] = npos[i, s, t]
    out["perm_pos_mean_null_a"], out["perm_pos_sd_null_a"], out["perm_pos_skew_null_a"] = pm[i, s, t], psd[i, s, t], pskew[i, s, t]
    if keep_top:
        flat = (i * K + s) * K + t
        where = np.searchsorted(scored, flat)
        hit = (where < len(scored)) & (scored[np.minimum(where, len(scored) - 1)] == flat)
        rows = np.full(len(out), None, dtype=object)
        for r, v in zip(np.flatnonzero(hit), -np.sort(-top[where[hit]], axis=1)):
            rows[r] = v
        out["perm_top_null_a"] = rows
    out.attrs["null_a"] = dict(version=version, B=int(B), seed=int(seed), keep_top=int(keep_top), rng=rng if permutations is None or rng == "cellchat" else "given",
                               data_max=data_max, n_lr=len(score.names), groups=[str(x) for x in levels])
    return out
