"""CellChat 2.2.0.9001 (jinworks/CellChat@75253cd), ported to Python. No R needed.

This is a line-by-line port of the parts of CellChat that produce ``net$prob`` and ``net$pval``
for single-cell RNA data, written to follow the R code as literally as possible. It is the
reference that ``nulls.null_a`` (a fast, vectorised Null A) is checked against.

    CellChatDB           read from the package's own data/CellChatDB.human.rda
    subsetDB             default: everything except "Non-protein Signaling"
    subsetData           data.signaling = the CellChatDB genes of the data
    identifyOverExpressedGenes          presto::wilcoxauc (b5df6ee) one-vs-rest, p < 0.05, logFC > 0
    identifyOverExpressedInteractions   variable.both = TRUE
    computeCommunProb    type = "triMean", Kh = 0.5, n = 1, nboot = 100, seed.use = 1,
                         population.size = FALSE
    set.seed / sample.int                R's Mersenne-Twister and rejection sampling (R >= 3.6), so
                                         the permutations are the ones CellChat draws

Sources: R/modeling.R, R/utilities.R, R/database.R of CellChat at 75253cd; R/wilcoxauc.R and
R/utils.R of presto at b5df6ee; src/main/RNG.c and random.c of R.
"""
from __future__ import annotations

import math
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import ndtr
from scipy.stats import rankdata

from ccc_nulls import datasets

COMMIT = "75253cd"
VERSION = "2.2.0.9001"
DB_URL = "https://raw.githubusercontent.com/jinworks/CellChat/{commit}/data/CellChatDB.{species}.rda"
ANNOTATION_LEVELS = ["Secreted Signaling", "ECM-Receptor", "Non-protein Signaling", "Cell-Cell Contact"]


# ====================================================================== R's random numbers
class RRandom:
    """R's default generator: ``set.seed(seed)`` then Mersenne-Twister ``unif_rand()``, and
    ``sample.int`` with ``sample.kind = "Rejection"`` (the default since R 3.6.0)."""

    I2_32M1 = 2.328306437080797e-10                    # 1 / (2^32 - 1)

    def __init__(self, seed: int):
        s = int(seed) & 0xFFFFFFFF
        for _ in range(50):                             # RNG_Init: initial scrambling
            s = (69069 * s + 1) & 0xFFFFFFFF
        key = np.empty(625, dtype=np.uint32)
        for j in range(625):                            # dummy[0] = mti, then mt[0..623]
            s = (69069 * s + 1) & 0xFFFFFFFF
            key[j] = s
        self._bg = np.random.MT19937()
        self._bg.state = {"bit_generator": "MT19937", "state": {"key": key[1:], "pos": 624}}  # FixupSeeds: mti = N
        self._buf = np.empty(0)
        self._i = 0

    def _refill(self, n=1 << 20):
        u = self._bg.random_raw(n).astype(np.float64) * 2.3283064365386963e-10   # MT_genrand
        u[u <= 0.0] = 0.5 * self.I2_32M1                                        # fixup()
        u[1.0 - u <= 0.0] = 1.0 - 0.5 * self.I2_32M1
        self._buf, self._i = u, 0

    def unif_rand(self) -> float:
        if self._i >= len(self._buf):
            self._refill()
        u = self._buf[self._i]
        self._i += 1
        return float(u)

    def runif(self, n: int) -> np.ndarray:
        return np.array([self.unif_rand() for _ in range(n)])

    def unif_index(self, dn: int) -> int:
        """R_unif_index: rejection sampling on the next power of two."""
        if dn <= 0:
            return 0
        bits = int(math.ceil(math.log2(dn)))
        mask = (1 << bits) - 1
        while True:
            v = 0
            for _ in range(0, bits + 1, 16):            # rbits()
                v = 65536 * v + int(math.floor(self.unif_rand() * 65536))
            dv = v & mask
            if dn > dv:
                return dv

    def sample_int(self, n: int, size: int | None = None) -> np.ndarray:
        """``sample.int(n, size)`` without replacement (do_sample, non-hash branch); 1-based."""
        size = n if size is None else size
        if n > 1e7 and size <= n / 2:
            raise NotImplementedError("R uses its hash sampler here")
        x = list(range(n))
        y = np.empty(size, dtype=np.int64)
        for i in range(size):
            j = self.unif_index(n)
            y[i] = x[j] + 1
            n -= 1
            x[j] = x[n]
        return y


def r_permutations(n_cells: int, nboot: int, seed: int) -> np.ndarray:
    """``set.seed(seed); replicate(nboot, sample.int(nC, size = nC))``: n_cells x nboot, 1-based."""
    rng = RRandom(seed)
    return np.column_stack([rng.sample_int(n_cells) for _ in range(nboot)])


# ====================================================================== CellChatDB
def db_dir(species: str = "human") -> Path:
    return datasets.data_dir() / "cellchat" / f"db_{species}_{COMMIT}"


def load_db(species: str = "human") -> dict:
    """CellChatDB.<species> from the package's data folder at the pinned commit, as DataFrames
    (cached as CSV under data/cellchat/). Needs the ``rdata`` package the first time."""
    d = db_dir(species)
    names = ["interaction", "complex", "cofactor", "geneInfo"]
    if not all((d / f"{n}.csv").exists() for n in names):
        import rdata

        d.mkdir(parents=True, exist_ok=True)
        rda = d / f"CellChatDB.{species}.rda"
        if not rda.exists():
            urllib.request.urlretrieve(DB_URL.format(commit=COMMIT, species=species), rda)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            obj = rdata.read_rda(str(rda))[f"CellChatDB.{species}"]
        for n in names:
            df = obj[n]
            df.index = df.index.astype(str)
            df.to_csv(d / f"{n}.csv")
    return {n: pd.read_csv(d / f"{n}.csv", index_col=0, keep_default_na=False, dtype=str) for n in names}


def subset_db(db: dict, non_protein: bool = False) -> dict:
    """``subsetDB(CellChatDB)`` with its defaults, then the annotation ordering of ``subsetData``."""
    search = ["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"] + (["Non-protein Signaling"] if non_protein else [])
    inter = db["interaction"]
    inter = inter[inter.annotation.isin(search)]
    if inter.annotation.nunique() > 1:                  # subsetData: order by annotation (stable)
        rank = inter.annotation.map({a: i for i, a in enumerate(ANNOTATION_LEVELS)})
        inter = inter.iloc[np.argsort(rank.to_numpy(), kind="stable")]
    return {**db, "interaction": inter}


def _subunits(complex_input: pd.DataFrame, name: str) -> list[str]:
    return [g for g in complex_input.loc[name].filter(like="subunit") if g]


def _cofactors(cofactor_input: pd.DataFrame, name: str) -> list[str]:
    return [g for g in cofactor_input.loc[name].filter(like="cofactor") if g]


def _extract_gene_subset(gene_set, complex_input, gene_info) -> list[str]:
    symbols = set(gene_info.Symbol)
    complexes = [g for g in gene_set if g not in symbols]
    genes = list(dict.fromkeys(g for g in gene_set if g in symbols))
    sub = complex_input.loc[[c for c in complexes if c in complex_input.index]].filter(like="subunit")
    subunits = [g for g in sub.to_numpy().ravel(order="F") if g]      # unlist(): column-major
    return list(dict.fromkeys(genes + subunits))


def extract_gene(db: dict) -> list[str]:
    """``extractGene``: every ligand, receptor, complex subunit and cofactor gene."""
    inter, cpx, cof, info = db["interaction"], db["complex"], db["cofactor"], db["geneInfo"]
    geneL = _extract_gene_subset(list(dict.fromkeys(inter.ligand)), cpx, info)
    geneR = _extract_gene_subset(list(dict.fromkeys(inter.receptor)), cpx, info)
    cofactor = list(dict.fromkeys(c for c in pd.concat([inter.agonist, inter.antagonist, inter.co_A_receptor,
                                                         inter.co_I_receptor]) if c))
    sub = cof.loc[[c for c in cofactor if c in cof.index]].filter(like="cofactor")
    gene_cof = [g for g in sub.to_numpy().ravel(order="F") if g]
    return list(dict.fromkeys(geneL + geneR + gene_cof))


# ====================================================================== presto::wilcoxauc
def wilcoxauc(X: sparse.spmatrix, codes: np.ndarray, K: int) -> dict:
    """presto's one-vs-rest Wilcoxon rank-sum test (normal approximation with continuity and
    tie correction) and its auxiliary statistics, for a cells x genes matrix of values >= 0."""
    X = sparse.csc_matrix(X, dtype=np.float64)
    N, G = X.shape
    n_g = np.bincount(codes, minlength=K).astype(np.float64)
    grs = np.empty((K, G))
    sums = np.empty((K, G))
    nnz = np.empty((K, G))
    tie_term = np.empty(G)
    for j in range(G):
        lo, hi = X.indptr[j], X.indptr[j + 1]
        rows, vals = X.indices[lo:hi], X.data[lo:hi]
        keep = vals != 0
        rows, vals = rows[keep], vals[keep]
        z = N - len(vals)
        r = rankdata(vals, method="average") + z if len(vals) else vals
        g = codes[rows]
        nz_g = np.bincount(g, minlength=K)
        grs[:, j] = (n_g - nz_g) * ((z + 1) / 2) + np.bincount(g, weights=r, minlength=K)
        sums[:, j] = np.bincount(g, weights=vals, minlength=K)
        nnz[:, j] = nz_g
        _, t = np.unique(vals, return_counts=True)
        t = np.append(t, z).astype(np.float64)
        tie_term[j] = np.sum(t ** 3 - t)
    ustat = grs - (n_g * (n_g + 1) / 2)[:, None]
    n1n2 = n_g * (N - n_g)
    zst = ustat - 0.5 * n1n2[:, None]
    zst = zst - np.sign(zst) * 0.5
    x1, x2 = float(N) ** 3 - N, 1 / (12 * (float(N) ** 2 - N))
    rhs = (x1 - tie_term) * x2
    usigma = np.sqrt(n1n2[:, None] * rhs[None, :])
    with np.errstate(invalid="ignore", divide="ignore"):
        zst = zst / usigma
    zst[~np.isfinite(zst)] = 0
    pval = 2 * ndtr(-np.abs(zst))
    means = sums / n_g[:, None]
    cs = sums.sum(0)
    lfc = means - (cs[None, :] - sums) / (N - n_g)[:, None]
    pct_in = 100 * nnz / n_g[:, None]
    pct_out = 100 * (nnz.sum(0)[None, :] - nnz) / (N - n_g)[:, None]
    return dict(pval=pval, logFC=lfc, pct_in=pct_in, pct_out=pct_out, avgExpr=means)


# ====================================================================== the pipeline
def labels_and_levels(labels) -> tuple[np.ndarray, np.ndarray]:
    labels = np.asarray(labels).astype(str)
    levels = np.unique(labels)
    return np.searchsorted(levels, labels), levels


def data_signaling(adata, db: dict):
    """``subsetData``: the CellChatDB genes of ``adata.X`` (cells x genes, float64 CSC), names."""
    genes = pd.Index(adata.var_names.astype(str))
    keep = np.where(genes.isin(set(extract_gene(db))))[0]
    X = sparse.csc_matrix(adata.X)[:, keep].astype(np.float64)
    if X.nnz and X.data.min() < 0:
        raise ValueError("CellChat expects non-negative normalised data")
    return X, list(genes[keep])


def identify_overexpressed_genes(X, codes, K, thresh_p=0.05, thresh_fc=0.0, thresh_pc=0.0) -> list[str] | np.ndarray:
    """``identifyOverExpressedGenes`` (do.fast = TRUE, only.pos = TRUE): gene indices significant
    and up in at least one group."""
    w = wilcoxauc(X, codes, K)
    pct_max = np.maximum(w["pct_in"], w["pct_out"])
    sig = (w["pval"] < thresh_p) & (np.abs(w["logFC"]) >= thresh_fc) & (pct_max > thresh_pc * 100) & (w["logFC"] > 0)
    return np.where(sig.any(0))[0]


def identify_overexpressed_interactions(db: dict, gene_use: list[str], features_sig: list[str]) -> pd.DataFrame:
    """``identifyOverExpressedInteractions`` (variable.both = TRUE): LRsig."""
    gene_use, features_sig = set(gene_use), set(features_sig)
    cpx = db["complex"]
    sig_complex = set()
    for name in cpx.index:
        sub = _subunits(cpx, name)
        if features_sig.intersection(sub) and all(g in gene_use for g in sub):
            sig_complex.add(name)
    ok = features_sig | sig_complex
    inter = db["interaction"]
    return inter[inter.ligand.isin(ok) & inter.receptor.isin(ok)].copy()


# ---------------------------------------------------------------------- computeCommunProb
def _r_mean(x: np.ndarray, axis=0) -> np.ndarray:
    """R's mean(): sequential sum / n, then one correction pass (long double == double on arm64)."""
    x = np.moveaxis(np.asarray(x, dtype=np.float64), axis, 0)
    n = x.shape[0]
    s = x[0].copy()
    for k in range(1, n):
        s = s + x[k]
    s = s / n
    with np.errstate(invalid="ignore"):
        t = np.zeros_like(s)
        for k in range(n):
            t = t + (x[k] - s)
        return np.where(np.isfinite(s), s + t / n, s)


def _r_prod(x: np.ndarray, axis=0) -> np.ndarray:
    x = np.moveaxis(np.asarray(x, dtype=np.float64), axis, 0)
    out = np.ones(x.shape[1:])
    for k in range(x.shape[0]):
        out = out * x[k]
    return out


def tri_mean_by_group(data_use: np.ndarray, codes: np.ndarray, K: int, quantile: str = "collapse") -> np.ndarray:
    """``aggregate(t(data.use), list(group), FUN = triMean)``: genes x groups.

    ``triMean`` is the mean of the quartiles (Q1, Q2, Q2, Q3), type 7:
    - ``quantile="collapse"`` (CellChat >= 2): ``collapse::fquantile`` interpolates ``a + h (b - a)``
      (fnth_fmedian_fquantile.c);
    - ``quantile="stats"`` (CellChat 1.x): ``stats::quantile`` interpolates ``(1 - h) a + h b``, and
      only where ``b != a`` (quantile.default).
    """
    data_use = np.atleast_2d(data_use)
    G = data_use.shape[0]
    out = np.empty((G, K))
    for k in range(K):
        x = np.sort(data_use[:, codes == k], axis=1)     # genes x n_k, each row sorted
        n = x.shape[1]
        q = []
        for Q in (0.25, 0.5, 0.75):
            h = (n - 1) * Q
            ih = int(h)
            frac = h - ih
            a = x[:, ih]
            if ih == n - 1 or frac <= 0.0:
                q.append(a)
            elif quantile == "collapse":
                q.append(a + frac * (x[:, ih + 1] - a))
            else:
                b = x[:, ih + 1]
                q.append(np.where(b != a, (1 - frac) * a + frac * b, a))
        out[:, k] = _r_mean(np.stack([q[0], q[1], q[1], q[2]]))
    return out


def geometric_mean(x: np.ndarray) -> np.ndarray:
    """``exp(apply(log(x), 2, mean))`` for a subunits x groups matrix."""
    with np.errstate(divide="ignore"):
        return np.exp(_r_mean(np.log(x), axis=0))


def compute_expr_lr(gene_lr: str, avg: pd.DataFrame, complex_input: pd.DataFrame) -> np.ndarray:
    """``computeExpr_LR`` for one ligand or receptor: its row of ``avg``, or for a complex the
    geometric mean of its subunits' rows (``computeExpr_complex``)."""
    if gene_lr in avg.index:
        return avg.loc[gene_lr].to_numpy()
    return geometric_mean(avg.loc[_subunits(complex_input, gene_lr)].to_numpy())


def compute_expr_coreceptor(cofactor_input, avg: pd.DataFrame, coreceptor: str) -> np.ndarray | float:
    """``computeExpr_coreceptor`` for one L-R pair."""
    if not coreceptor:
        return np.ones(avg.shape[1])
    genes = [g for g in _cofactors(cofactor_input, coreceptor) if g in avg.index]
    if len(genes) == 1:
        return 1 + avg.loc[genes[0]].to_numpy()
    if len(genes) > 1:
        return _r_prod(1 + avg.loc[genes].to_numpy(), axis=0)
    return np.ones(avg.shape[1])


def compute_expr_agonist(avg: pd.DataFrame, cofactor_input, agonist: str, Kh, n) -> np.ndarray:
    genes = [g for g in _cofactors(cofactor_input, agonist) if g in avg.index]
    a = avg.loc[genes].to_numpy() if genes else None
    if len(genes) == 1:
        return (1 + a ** n / (Kh ** n + a ** n))[0] if n != 1 else (1 + a / (Kh + a))[0]
    if len(genes) > 1:
        return _r_prod(1 + a / (Kh + a) if n == 1 else 1 + a ** n / (Kh ** n + a ** n), axis=0)
    return np.ones(avg.shape[1])


def compute_expr_antagonist(avg: pd.DataFrame, cofactor_input, antagonist: str, Kh, n) -> np.ndarray:
    genes = [g for g in _cofactors(cofactor_input, antagonist) if g in avg.index]
    a = avg.loc[genes].to_numpy() if genes else None
    if len(genes) == 1:
        return (Kh / (Kh + a))[0] if n == 1 else (Kh ** n / (Kh ** n + a ** n))[0]
    if len(genes) > 1:
        return _r_prod(Kh / (Kh + a) if n == 1 else Kh ** n / (Kh ** n + a ** n), axis=0)
    return np.ones(avg.shape[1])


def compute_commun_prob(X, genes, codes, levels, LRsig: pd.DataFrame, db: dict, *, nboot=100, seed=1,
                        Kh=0.5, n=1, population_size=False, permutation=None, log=None,
                        version: str = "2.2") -> dict:
    """``computeCommunProb(type = "triMean")`` for scRNA-seq data, following the R loop.

    ``X`` is data.signaling (cells x signaling genes). Returns ``prob`` and ``pval``
    (nLR x K x K), the permutation matrix and the scaling constant.
    ``version="2.2"``: CellChat 2.2.0.9001 (R/modeling.R at 75253cd). ``version="1.0"``: CellChat
    1.0.0, the version of the 2021 paper (sqjin/CellChat@v1.0.0), run on the raw data
    (``raw.use = TRUE``; its default projected the data on a protein network first)."""
    if version == "1.0":
        return _compute_commun_prob_v1(X, genes, codes, levels, LRsig, db, nboot=nboot, seed=seed, Kh=Kh,
                                       n=n, population_size=population_size, permutation=permutation, log=log)
    if version != "2.2":
        raise ValueError("version must be '2.2' or '1.0'")
    complex_input, cofactor_input = db["complex"], db["cofactor"]
    K, nC = len(levels), X.shape[0]
    data_max = float(X.max())
    # only the genes LRsig touches are needed; data.use = data / max(data) over all signaling genes
    need, gene_set = set(), set(genes)
    for _, r in LRsig.iterrows():
        for x in (r.ligand, r.receptor):
            need.update([x] if x in gene_set else _subunits(complex_input, x))
        for c in (r.agonist, r.antagonist, r.co_A_receptor, r.co_I_receptor):
            if c:
                need.update(_cofactors(cofactor_input, c))
    gpos = {g: j for j, g in enumerate(genes)}
    need = [g for g in genes if g in need]
    data_use = (X[:, [gpos[g] for g in need]].toarray() / data_max).T          # genes x cells

    def avg_for(cd):
        return pd.DataFrame(tri_mean_by_group(data_use, cd, K), index=need)

    avg = avg_for(codes)
    if permutation is None:
        permutation = r_permutations(nC, nboot, seed)
    avg_boot = [avg_for(codes[permutation[:, b] - 1]) for b in range(nboot)]
    if log:
        log("  group averages for", nboot, "permutations done")
    sizes = np.bincount(codes, minlength=K) / nC

    def p1_of(r, av):
        L = compute_expr_lr(r.ligand, av, complex_input)
        R = compute_expr_lr(r.receptor, av, complex_input)
        R = R * compute_expr_coreceptor(cofactor_input, av, r.co_A_receptor) / \
            compute_expr_coreceptor(cofactor_input, av, r.co_I_receptor)
        LR = np.outer(L, R)
        return LR / (Kh + LR) if n == 1 else LR ** n / (Kh ** n + LR ** n)

    nLR = len(LRsig)
    Prob = np.zeros((nLR, K, K))
    Pval = np.zeros((nLR, K, K))
    for i, (_, r) in enumerate(LRsig.iterrows()):
        P1 = p1_of(r, avg)
        if P1.sum() == 0:
            Prob[i] = P1
            Pval[i] = 1
            continue

        def extra(av):
            P2 = np.outer(*(2 * [compute_expr_agonist(av, cofactor_input, r.agonist, Kh, n)])) if r.agonist else 1.0
            P3 = np.outer(*(2 * [compute_expr_antagonist(av, cofactor_input, r.antagonist, Kh, n)])) if r.antagonist else 1.0
            return P2, P3

        P2, P3 = extra(avg)
        P4 = np.outer(sizes, sizes) if population_size else 1.0
        Pnull = P1 * P2 * P3 * P4
        Prob[i] = Pnull
        n_reject = np.zeros((K, K), dtype=np.int64)
        for b in range(nboot):
            P2b, P3b = extra(avg_boot[b])
            Pboot = p1_of(r, avg_boot[b]) * P2b * P3b * P4
            n_reject += (Pboot - Pnull) > 0
        Pval[i] = n_reject / nboot
    Pval[Prob == 0] = 1
    return dict(prob=Prob, pval=Pval, permutation=permutation, data_max=data_max, groups=list(levels),
                LRsig=LRsig, genes_signaling=genes)


# ---------------------------------------------------------------------- computeCommunProb, CellChat 1.0.0
def cell_expr_lr(gene_lr: str, data_use: pd.DataFrame, complex_input: pd.DataFrame) -> np.ndarray:
    """v1.0.0 ``computeExpr_LR`` / ``computeExpr_complex`` per cell: the gene's row of ``data.use``
    (genes x cells), or for a complex ``geometricMean(data.use[subunits, ])`` across subunits."""
    if gene_lr in data_use.index:
        return data_use.loc[gene_lr].to_numpy()
    return geometric_mean(data_use.loc[_subunits(complex_input, gene_lr)].to_numpy())


def _compute_commun_prob_v1(X, genes, codes, levels, LRsig, db, *, nboot, seed, Kh, n, population_size,
                            permutation, log):
    """CellChat 1.0.0 ``computeCommunProb`` (R/modeling.R at sqjin/CellChat@v1.0.0), raw data."""
    if population_size or n != 1:
        raise NotImplementedError("version 1.0 is ported for population.size = FALSE and n = 1")
    complex_input, cofactor_input = db["complex"], db["cofactor"]
    K, nC = len(levels), X.shape[0]
    data_max = float(X.max())
    need, gene_set = set(), set(genes)
    for _, r in LRsig.iterrows():
        for x in (r.ligand, r.receptor):
            need.update([x] if x in gene_set else _subunits(complex_input, x))
        for c in (r.agonist, r.antagonist, r.co_A_receptor, r.co_I_receptor):
            if c:
                need.update(_cofactors(cofactor_input, c))
    gpos = {g: j for j, g in enumerate(genes)}
    need = [g for g in genes if g in need]
    data_use = pd.DataFrame((X[:, [gpos[g] for g in need]].toarray() / data_max).T, index=need)  # genes x cells
    if permutation is None:
        permutation = r_permutations(nC, nboot, seed)
    boot_codes = [codes[permutation[:, b] - 1] for b in range(nboot)]

    # agonist / antagonist genes: per-gene trimean by group (computeExpr_agonist / _antagonist)
    cof_genes = sorted({g for _, r in LRsig.iterrows() for c in (r.agonist, r.antagonist) if c
                        for g in _cofactors(cofactor_input, c) if g in data_use.index})
    tm = lambda M, cd: tri_mean_by_group(M, cd, K, quantile="stats")
    cof = data_use.loc[cof_genes].to_numpy()
    avg_of = lambda cd: (pd.DataFrame(tm(cof, cd), index=cof_genes) if cof_genes
                         else pd.DataFrame(np.zeros((0, K))))
    cof_avg = avg_of(codes)
    cof_boot = [avg_of(cd) for cd in boot_codes]

    def extra(r, av):
        P2 = np.outer(*(2 * [compute_expr_agonist(av, cofactor_input, r.agonist, Kh, n)])) if r.agonist else 1.0
        P3 = np.outer(*(2 * [compute_expr_antagonist(av, cofactor_input, r.antagonist, Kh, n)])) if r.antagonist else 1.0
        return P2, P3

    nLR = len(LRsig)
    Prob = np.zeros((nLR, K, K))
    Pval = np.zeros((nLR, K, K))
    for i, (_, r) in enumerate(LRsig.iterrows()):
        dataL = cell_expr_lr(r.ligand, data_use, complex_input)
        dataR = cell_expr_lr(r.receptor, data_use, complex_input)
        dataR = dataR * compute_expr_coreceptor(cofactor_input, data_use, r.co_A_receptor) / \
            compute_expr_coreceptor(cofactor_input, data_use, r.co_I_receptor)
        LR = np.outer(tm(dataL, codes)[0], tm(dataR, codes)[0])
        P1 = LR / (Kh + LR)
        if P1.sum() == 0:
            Prob[i] = P1
            Pval[i] = 1
            continue
        P2, P3 = extra(r, cof_avg)
        Pnull = P1 * P2 * P3
        Prob[i] = Pnull
        n_reject = np.zeros((K, K), dtype=np.int64)
        for b in range(nboot):
            LRb = np.outer(tm(dataL, boot_codes[b])[0], tm(dataR, boot_codes[b])[0])
            P2b, P3b = extra(r, cof_boot[b])
            n_reject += (LRb / (Kh + LRb) * P2b * P3b - Pnull) >= 0       # v1.0.0: >=
        Pval[i] = n_reject / nboot
        if log and (i + 1) % 500 == 0:
            log(f"  v1.0: {i + 1}/{nLR} L-R pairs")
    Pval[Prob == 0] = 1
    return dict(prob=Prob, pval=Pval, permutation=permutation, data_max=data_max, groups=list(levels),
                LRsig=LRsig, genes_signaling=genes, version="1.0")
