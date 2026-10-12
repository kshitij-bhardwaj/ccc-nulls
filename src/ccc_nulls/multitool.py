"""Null-model comparison across LIANA-py's scoring methods.

Every LIANA score is a fixed function of a few per-cell-type statistics of the ligand
(in the source type) and the receptor (in the target type):

    mean, fraction expressing, z-score, one-vs-rest log2FC, trimean, sum over types.

:class:`StatsEngine` recomputes those statistics for any labelling of the cells (real,
permuted, permuted within patient). :func:`tool_scores` applies each tool's formula, mirroring
liana-py 1.x line for line, and is validated against LIANA's own output in
``benchmarks/tnbc_openproblems.py``. The nulls then differ only in *what* is randomised:
cell labels (N1, N2, N5) or gene identity (N3, N4).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats as sstats

# tool -> (liana-py method name, statistic tested (LIANA column), statistics it needs)
TOOLS = {
    "CellPhoneDB": ("cellphonedb", "lr_means", {"means"}),
    "GeometricMean": ("geometric_mean", "lr_gmeans", {"means"}),
    "CellChat": ("cellchat", "lr_probs", {"trimean"}),
    "Connectome": ("connectome", "scaled_weight", {"zscores"}),
    "NATMI": ("natmi", "spec_weight", {"means", "sums"}),
    "SingleCellSignalR": ("singlecellsignalr", "lrscore", {"means"}),
    "log2FC": ("logfc", "lr_logfc", {"logfc"}),
}
CELLCHAT_KH = 0.5


# ------------------------------------------------------------------ statistics
def _indicator(codes: np.ndarray, K: int) -> sparse.csr_matrix:
    n = len(codes)
    return sparse.csr_matrix((np.ones(n), (codes, np.arange(n))), shape=(K, n))


def _quantile_sorted_with_zeros(d_sorted, n, q):
    """numpy-style linear quantile of a vector of n values = (n - len(d)) zeros + d_sorted."""
    nz0 = n - len(d_sorted)
    h = (n - 1) * q
    lo = int(np.floor(h))
    hi = min(lo + 1, n - 1)
    v_lo = 0.0 if lo < nz0 else d_sorted[lo - nz0]
    v_hi = 0.0 if hi < nz0 else d_sorted[hi - nz0]
    return v_lo + (h - lo) * (v_hi - v_lo)


def trimean_by_type(Xcsr: sparse.csr_matrix, codes: np.ndarray, K: int, mat_max: float) -> np.ndarray:
    """K x G trimean of X / mat_max per cell type, exactly as liana's _trimean
    ((Q1 + 2 * median + Q3) / 4, numpy linear quantiles), using sparsity."""
    G = Xcsr.shape[1]
    out = np.zeros((K, G))
    for k in range(K):
        Xk = Xcsr[codes == k].tocsc()
        n = Xk.shape[0]
        if n == 0:
            continue
        ip, data = Xk.indptr, Xk.data
        nnz = np.diff(ip)
        # if fewer than ~25% of cells express the gene, Q1 = median = Q3 = 0
        cand = np.where(nnz > 0)[0]
        h3 = (n - 1) * 0.75
        need = cand[(int(np.floor(h3)) + 1) >= (n - nnz[cand])]  # Q3 can touch non-zeros
        for j in need:
            d = np.sort(data[ip[j]:ip[j + 1]]) / mat_max
            q1 = _quantile_sorted_with_zeros(d, n, 0.25)
            q2 = _quantile_sorted_with_zeros(d, n, 0.5)
            q3 = _quantile_sorted_with_zeros(d, n, 0.75)
            out[k, j] = (q1 + 2 * q2 + q3) / 4
    return out


class StatsEngine:
    """Per-cell-type statistics for a fixed set of gene columns, under any labelling.

    X is the normalised (log1p, natural base) cells x genes matrix LIANA scored.
    mat_mean / mat_max are LIANA's global constants (whole matrix), which no null changes.
    Scaling for z-scores uses all cells (as sc.pp.scale does), so it is also fixed.
    """

    def __init__(self, X, genes_idx: np.ndarray, mat_mean: float, mat_max: float, scale_ddof: int = 1):
        # float32 storage keeps whole-transcriptome engines within laptop memory; statistics
        # are accumulated in float64 by the sparse products below.
        Xs = sparse.csr_matrix(X)[:, genes_idx].astype(np.float32)
        self.X = Xs.tocsr()
        self.Xbin = (self.X > 0).astype(np.float32).tocsr()
        E = self.X.copy()
        E.data = np.expm1(E.data)
        self.E = E.tocsr()
        n = self.X.shape[0]
        mu = np.asarray(self.X.mean(0, dtype=np.float64)).ravel()
        sq = np.asarray(self.X.multiply(self.X).mean(0, dtype=np.float64)).ravel()
        var = np.maximum(sq - mu**2, 0) * (n / (n - scale_ddof) if scale_ddof else 1.0)
        sd = np.sqrt(var)
        sd[sd == 0] = 1.0
        self.mu, self.sd = mu, sd
        self.E_tot = np.asarray(self.E.sum(0, dtype=np.float64)).ravel()
        self.n = n
        self.mat_mean, self.mat_max = mat_mean, mat_max

    def compute(self, codes: np.ndarray, K: int, need: set) -> dict:
        O = _indicator(codes, K)  # float64, so O @ X accumulates in float64
        size = np.asarray(O.sum(1)).ravel()[:, None]
        S = {}
        if need & {"means", "zscores", "sums", "props"}:
            means = (O @ self.X).toarray() / size
            S["means"] = means
            S["zscores"] = (means - self.mu) / self.sd
            S["sums"] = means.sum(0, keepdims=True)
        if "props" in need:
            S["props"] = (O @ self.Xbin).toarray() / size
        if "logfc" in need:
            Ek = (O @ self.E).toarray()
            subj = Ek / size
            rest = (self.E_tot - Ek) / np.maximum(self.n - size, 1)
            S["logfc"] = np.log2(subj + 1) - np.log2(rest + 1)
        if "trimean" in need:
            S["trimean"] = trimean_by_type(self.X, codes, K, self.mat_max)
        return S


def tool_scores(tool: str, S: dict, li, ri, si, ti, mat_mean: float) -> np.ndarray:
    """Each tool's statistic, mirroring liana-py's method functions."""
    if tool == "CellPhoneDB":
        L, R = S["means"][si, li], S["means"][ti, ri]
        out = (L + R) / 2
        out[(L == 0) | (R == 0)] = 0
        return out
    if tool == "GeometricMean":
        L, R = S["means"][si, li], S["means"][ti, ri]
        return np.sqrt(L * R)
    if tool == "CellChat":
        p = S["trimean"][si, li] * S["trimean"][ti, ri]
        return p / (CELLCHAT_KH + p)
    if tool == "Connectome":
        return (S["zscores"][si, li] + S["zscores"][ti, ri]) / 2
    if tool == "NATMI":
        m, s = S["means"], S["sums"][0]
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.nan_to_num((m[si, li] / s[li]) * (m[ti, ri] / s[ri]))
    if tool == "SingleCellSignalR":
        lr = np.sqrt(S["means"][si, li]) * np.sqrt(S["means"][ti, ri])
        return lr / (lr + mat_mean)
    if tool == "log2FC":
        return (S["logfc"][si, li] + S["logfc"][ti, ri]) / 2
    raise ValueError(tool)


# ------------------------------------------------------------------ candidates
class Candidates:
    """One tool's candidate rows, as index arrays into an engine's gene space."""

    def __init__(self, df: pd.DataFrame, gene_pos: dict, type_pos: dict):
        self.df = df.reset_index(drop=True)
        self.li = df.ligand.map(gene_pos).to_numpy()
        self.ri = df.receptor.map(gene_pos).to_numpy()
        self.si = df.source.map(type_pos).to_numpy()
        self.ti = df.target.map(type_pos).to_numpy()

    def __len__(self):
        return len(self.df)


# ------------------------------------------------------------------ label nulls
def permute_codes(codes, rng, groups=None):
    """Uniform permutation (groups=None) or permutation within each group (e.g. patient)."""
    if groups is None:
        return rng.permutation(codes)
    out = codes.copy()
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        out[idx] = codes[rng.permutation(idx)]
    return out


def label_null(engine: StatsEngine, codes, K, cands: dict, obs: dict, B: int, rng, groups=None,
               tools=None, log_every=0) -> dict:
    """N1 (groups=None) / N2 (groups=patient). Returns p-values per tool."""
    tools = tools or list(cands)
    need = set().union(*(TOOLS[t][2] for t in tools))
    exceed = {t: np.zeros(len(cands[t])) for t in tools}
    for b in range(B):
        S = engine.compute(permute_codes(codes, rng, groups), K, need)
        for t in tools:
            c = cands[t]
            s = tool_scores(t, S, c.li, c.ri, c.si, c.ti, engine.mat_mean)
            exceed[t] += s >= obs[t] - 1e-12
        if log_every and (b + 1) % log_every == 0:
            print(f"    perm {b + 1}/{B}", flush=True)
    return {t: (1 + exceed[t]) / (B + 1) for t in tools}


def patient_replicate_null(X, codes, K, patients, gene_idx, mat_mean, mat_max, cands: dict, B: int, rng,
                           min_cells: int = 10, min_patients: int = 3, tools=None) -> dict:
    """N5: per-patient z-scores (within-patient permutation), one-sided t-test across patients."""
    tools = tools or list(cands)
    need = set().union(*(TOOLS[t][2] for t in tools))
    Z = {t: [] for t in tools}
    Xc = sparse.csr_matrix(X)
    for p in np.unique(patients):
        cells = np.where(patients == p)[0]
        eng = StatsEngine(Xc[cells], gene_idx, mat_mean, mat_max)
        cp = codes[cells]
        present = np.bincount(cp, minlength=K) >= min_cells
        S = eng.compute(cp, K, need)
        perm = {t: [] for t in tools}
        for _ in range(B):
            Sb = eng.compute(rng.permutation(cp), K, need)
            for t in tools:
                c = cands[t]
                perm[t].append(tool_scores(t, Sb, c.li, c.ri, c.si, c.ti, mat_mean))
        for t in tools:
            c = cands[t]
            o = tool_scores(t, S, c.li, c.ri, c.si, c.ti, mat_mean)
            P = np.vstack(perm[t])
            sd = P.std(0)
            z = (o - P.mean(0)) / np.where(sd > 0, sd, np.nan)
            z[~(present[c.si] & present[c.ti])] = np.nan
            Z[t].append(z)
    out = {}
    for t in tools:
        Zt = np.vstack(Z[t]).T                        # rows x patients
        n = np.sum(~np.isnan(Zt), 1)
        m = np.nanmean(Zt, 1)
        s = np.nanstd(Zt, 1, ddof=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            tstat = m / (s / np.sqrt(n))
        p = sstats.t.sf(tstat, df=np.maximum(n - 1, 1))
        p[(n < min_patients) | ~np.isfinite(p)] = 1.0
        out[t] = p
    return out


# ------------------------------------------------------------------ gene-identity nulls
def decoy_null(S_all: dict, mat_mean, cands_all: dict, obs: dict, strata: np.ndarray, pool_mask: np.ndarray,
               props: np.ndarray, expr_prop: float, B: int, rng, matched: bool = True, tools=None,
               chunk: int | None = None) -> dict:
    """N3 (matched=True) / N4 (matched=False). Monte-Carlo over the reference set
    (decoys + the candidate's own gene); p = (1 + #{draws >= obs}) / (B + 1) is valid because the
    candidate is itself a uniform member of the set it is compared with.

    S_all : statistics for ALL genes under the real labels; cands_all index into all genes.
    props : K x G fraction expressing (real labels), for the expression filter on decoys.
    """
    tools = tools or list(cands_all)
    K, G = props.shape
    chunk = chunk or max(50, 4_000_000 // B)        # keep each chunk at ~4M draws
    pools = {}

    def pool(stratum, k):
        key = (stratum if matched else -1, k)
        if key not in pools:
            ok = pool_mask & (props[k] >= expr_prop)
            if matched:
                ok &= strata == stratum
            pools[key] = np.where(ok)[0]
        return pools[key]

    out = {}
    for t in tools:
        c = cands_all[t]
        n = len(c)
        exceed = np.zeros(n)
        for a in range(0, n, chunk):
            sl = slice(a, min(a + chunk, n))
            li, ri, si, ti = c.li[sl], c.ri[sl], c.si[sl], c.ti[sl]
            m = len(li)
            Ld = np.empty((m, B), dtype=np.int64)
            Rd = np.empty((m, B), dtype=np.int64)
            for j in range(m):
                pl = pool(strata[li[j]], si[j])
                pr = pool(strata[ri[j]], ti[j])
                il = rng.integers(0, len(pl) + 1, B)          # index len(pl) = the ligand itself
                ir = rng.integers(0, len(pr) + 1, B)
                Ld[j] = np.where(il == len(pl), li[j], pl[np.minimum(il, max(len(pl) - 1, 0))] if len(pl) else li[j])
                Rd[j] = np.where(ir == len(pr), ri[j], pr[np.minimum(ir, max(len(pr) - 1, 0))] if len(pr) else ri[j])
            s = tool_scores(t, S_all, Ld.ravel(), Rd.ravel(), np.repeat(si, B), np.repeat(ti, B), mat_mean).reshape(m, B)
            exceed[sl] = (s >= obs[t][sl, None] - 1e-12).sum(1)
        out[t] = (1 + exceed) / (B + 1)
    return out
