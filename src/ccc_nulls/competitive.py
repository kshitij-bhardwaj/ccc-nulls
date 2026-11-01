"""Null B2: the decoy-gene (competitive) null with the inter-gene correlation put back into its variance.

The October decoy null (:func:`ccc_nulls.nulls.null_b`) ranked a candidate against abundance-matched
genes outside the ligand-receptor databases, as if the decoys were independent draws from the
candidate's null. Genes are correlated across cells, so that premise is false, and the test was
anti-conservative (p < 0.05 for 13-14% of tests on shuffled labels). This is the known flaw of
competitive gene-set tests by gene sampling (Goeman & Buehlmann 2007, Bioinformatics 23:980). CAMERA
(Wu & Smyth 2012, Nucleic Acids Res. 40:e133) keeps the competitive null and puts the inter-gene
correlation into the variance; Efron & Tibshirani's restandardization (2007, Ann. Appl. Stat. 1:107)
is the same idea.

One side of an edge (the ligand's subunit genes S in the source type, or the receptor's in the
target type) is tested against its decoys D:

    t_g          one-vs-rest Welch t of gene g in the type (ln(1 + CP10k))
    T            mean_S t - mean_D t
    Var(T)       tau2 (VIF_S/m + 1/n) + VIF_S/m + VIF_D/n - 2 r_SD
    VIF_S        1 + (m - 1) rho_S            (CAMERA's variance inflation factor)
    VIF_D        1 + (n - 1) rho_D
    tau2         max(0, var_D(t) - (1 - rho_D))   spread of real specificity among comparable genes

rho_S, rho_D and r_SD are mean Pearson correlations of within-type residuals (CAMERA's residual
correlation) among S, among D, and between S and D. The side's p is the upper tail of t with n - 1
degrees of freedom at T / sqrt(Var(T)), and an edge's p is max(p_ligand, p_receptor), the
intersection-union rule of :mod:`ccc_nulls.composite`. ``correction=False`` gives the old premise,
Var(T) = var_D(t) (1/m + 1/n). Plan and validity checks: benchmarks/liana_nulls/AMENDMENT_3.md.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats as sstats

from .nulls import gene_strata

MIN_FRAC = 0.10
MAX_DECOYS = 200
N_CELLS_CORR = 5000


@dataclass
class GeneStats:
    """What every side test of one data set reads: per-type t and detection, strata, residuals."""

    genes: pd.Index
    types: list
    t: np.ndarray            # types x genes, one-vs-rest Welch t
    frac: np.ndarray         # types x genes, fraction of the type's cells with the gene detected
    strata: np.ndarray       # genes
    U: np.ndarray | None = None         # cells (subsample) x used genes, within-type residuals, unit norm
    u_col: np.ndarray | None = None     # gene -> column of U (-1: not computed)

    def with_residuals(self, X, codes, genes_used, n_cells: int = N_CELLS_CORR, seed: int = 0) -> "GeneStats":
        """Residuals for ``genes_used`` (indices) only: decoys and candidates, not every gene (memory)."""
        used = np.unique(np.asarray(genes_used, dtype=np.int64))
        self.u_col = np.full(len(self.genes), -1, dtype=np.int64)
        self.u_col[used] = np.arange(len(used))
        self.U = residuals(X, np.asarray(codes), len(self.types), n_cells, seed, cols=used)
        return self

    def u(self, idx) -> np.ndarray:
        return self.U[:, self.u_col[np.asarray(idx)]]


def welch_t(X, codes: np.ndarray, n_types: int) -> tuple[np.ndarray, np.ndarray]:
    """(t, frac): one-vs-rest Welch t and detection fraction, types x genes."""
    X = sparse.csr_matrix(X, dtype=np.float64) if sparse.issparse(X) else np.asarray(X, np.float64)
    n = np.bincount(codes, minlength=n_types).astype(float)
    ind = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(n_types, len(codes)))
    dense = lambda M: M.toarray() if sparse.issparse(M) else np.asarray(M)
    s1 = dense(ind @ X)
    s2 = dense(ind @ (X.multiply(X) if sparse.issparse(X) else X * X))
    nz = dense(ind @ ((X > 0).astype(np.float64)))
    N, S1, S2 = n.sum(), s1.sum(0), s2.sum(0)
    nr = N - n
    m_k, m_r = s1 / n[:, None], (S1 - s1) / nr[:, None]
    v_k = np.maximum(s2 - n[:, None] * m_k**2, 0) / np.maximum(n - 1, 1)[:, None]
    v_r = np.maximum((S2 - s2) - nr[:, None] * m_r**2, 0) / np.maximum(nr - 1, 1)[:, None]
    se = np.sqrt(v_k / n[:, None] + v_r / nr[:, None])
    t = np.divide(m_k - m_r, se, out=np.zeros_like(se), where=se > 0)
    return t, nz / n[:, None]


def residuals(X, codes: np.ndarray, n_types: int, n_cells: int = N_CELLS_CORR, seed: int = 0, cols=None) -> np.ndarray:
    """Within-type residuals of a cell subsample, each gene's column scaled to unit norm, so that
    ``U[:, g] @ U[:, h]`` is the residual correlation of g and h (0 for a constant gene)."""
    rng = np.random.default_rng(seed)
    cells = np.sort(rng.choice(len(codes), n_cells, replace=False)) if len(codes) > n_cells else np.arange(len(codes))
    Xs = X[cells]
    if cols is not None:
        Xs = Xs[:, cols]
    Xs = np.asarray(Xs.todense() if sparse.issparse(Xs) else Xs, dtype=np.float32)
    c = codes[cells]
    for k in np.unique(c):
        rows = c == k
        Xs[rows] -= Xs[rows].mean(0)
    norm = np.linalg.norm(Xs, axis=0)
    Xs /= np.where(norm > 0, norm, 1.0)
    return Xs


def gene_stats(X, genes, codes, n_types: int, types=None) -> GeneStats:
    """t, detection and strata for cells x genes ``X`` (ln(1 + CP10k)) and type codes; residuals are
    added by :meth:`GeneStats.with_residuals` once the decoys are known."""
    codes = np.asarray(codes)
    t, frac = welch_t(X, codes, n_types)
    return GeneStats(genes=pd.Index(np.asarray(genes, dtype=str)), types=list(types) if types is not None else list(range(n_types)),
                     t=t, frac=frac, strata=gene_strata(X))


def decoy_groups(gs: GeneStats, excluded, matched: bool = True, min_frac: float = MIN_FRAC,
                 max_decoys: int = MAX_DECOYS, seed: int = 0) -> dict:
    """Decoys per (stratum, type), or per type when ``matched`` is False (key (-1, type)): genes not in
    ``excluded`` (the databases), detected in >= min_frac of the type's cells, at most max_decoys."""
    rng = np.random.default_rng(seed)
    pool = ~pd.Index(gs.genes.str.upper()).isin({str(g).upper() for g in excluded})     # either species' case
    out = {}
    strata = np.unique(gs.strata) if matched else [-1]
    for st in strata:
        in_st = (gs.strata == st) if matched else np.ones(len(gs.genes), bool)
        for k in range(len(gs.types)):
            g = np.where(in_st & pool & (gs.frac[k] >= min_frac))[0]
            if len(g) > max_decoys:
                g = np.sort(rng.choice(g, max_decoys, replace=False))
            out[(int(st), k)] = g
    return out


def _mean_corr(u_sum: np.ndarray, diag: float, n: int) -> float:
    return float((u_sum @ u_sum - diag) / (n * (n - 1))) if n > 1 else 0.0


def side_tests(gs: GeneStats, sides: pd.DataFrame, groups: dict, matched: bool = True,
               min_frac: float = MIN_FRAC) -> pd.DataFrame:
    """One row per side (columns ``complex``: subunit symbols joined by "_", ``type``: code), with the
    gate, the CAMERA-corrected z and p, the uncorrected z and p, and the pieces of the variance."""
    gidx = {g: i for i, g in enumerate(gs.genes)}
    ref_cache = {}
    out = []
    for cx, k in zip(sides["complex"].to_numpy(), sides["type"].to_numpy()):
        sub = [gidx.get(s) for s in str(cx).split("_")]
        rec = dict(gate=False, m=len(sub), n=0, T=np.nan, z=np.nan, p=1.0, z_uncorrected=np.nan, p_uncorrected=1.0,
                   rho_S=np.nan, rho_D=np.nan, r_SD=np.nan, tau2=np.nan)
        if any(i is None for i in sub) or (gs.frac[k, sub] < min_frac).any():
            out.append(rec)
            continue
        key = (tuple(sorted({int(gs.strata[i]) for i in sub})) if matched else (-1,), int(k))
        if key not in ref_cache:
            D = np.unique(np.concatenate([groups[(st, int(k))] for st in key[0]]))
            UD = gs.u(D).astype(np.float64)
            uD, diag = UD.sum(1), float((UD**2).sum())
            tD = gs.t[k, D]
            ref_cache[key] = (D, uD, _mean_corr(uD, diag, len(D)), tD.mean() if len(D) else np.nan,
                              tD.var(ddof=1) if len(D) > 1 else np.nan)
        D, uD, rho_D, mD, vD = ref_cache[key]
        n, m = len(D), len(sub)
        rec["n"] = n
        if n < 3:
            out.append(rec)
            continue
        US = gs.u(sub).astype(np.float64)
        uS = US.sum(1)
        rho_S = _mean_corr(uS, float((US**2).sum()), m)
        r_SD = float(uS @ uD / (m * n))
        T = float(gs.t[k, sub].mean() - mD)
        tau2 = max(0.0, vD - (1 - rho_D))
        vif_S, vif_D = 1 + (m - 1) * rho_S, 1 + (n - 1) * rho_D
        var = tau2 * (vif_S / m + 1 / n) + vif_S / m + vif_D / n - 2 * r_SD
        var_u = vD * (1 / m + 1 / n)
        z = T / np.sqrt(var) if var > 0 else np.nan
        zu = T / np.sqrt(var_u) if var_u > 0 else np.nan
        rec.update(gate=True, T=T, z=z, p=float(sstats.t.sf(z, n - 1)) if np.isfinite(z) else 1.0,
                   z_uncorrected=zu, p_uncorrected=float(sstats.t.sf(zu, n - 1)) if np.isfinite(zu) else 1.0,
                   rho_S=rho_S, rho_D=rho_D, r_SD=r_SD, tau2=tau2)
        out.append(rec)
    return pd.concat([sides.reset_index(drop=True), pd.DataFrame(out)], axis=1)


VARIANTS = ("b2", "b2_uncorrected", "b2_unmatched", "b2_joint")


def null_b2(X, gs: GeneStats, codes, rows: pd.DataFrame, excluded, seed: int = 0,
            n_cells_corr: int = N_CELLS_CORR) -> pd.DataFrame:
    """Null B2 and its ablation variants for rows with ``ligand_complex, receptor_complex, source,
    target`` (type names as in ``gs.types``). Returns one row per input row: ``p_L_b2, p_R_b2,
    p_b2`` (full), ``p_b2_uncorrected`` (no correlation), ``p_b2_unmatched`` (any non-database
    decoy), ``p_b2_joint`` (z of both sides combined), and ``gate_b2``."""
    code = {str(t): i for i, t in enumerate(gs.types)}
    s = np.array([code.get(str(x), -1) for x in rows.source])
    t = np.array([code.get(str(x), -1) for x in rows.target])
    if (s < 0).any() or (t < 0).any():
        raise ValueError("rows name cell types that are not in the gene statistics")
    lig = pd.DataFrame({"complex": rows.ligand_complex.astype(str).to_numpy(), "type": s})
    rec = pd.DataFrame({"complex": rows.receptor_complex.astype(str).to_numpy(), "type": t})
    sides = pd.concat([lig, rec]).drop_duplicates().reset_index(drop=True)
    groups = decoy_groups(gs, excluded, matched=True, seed=seed)
    groups_u = decoy_groups(gs, excluded, matched=False, seed=seed)
    gidx = {g: i for i, g in enumerate(gs.genes)}
    cand = [gidx[g] for cx in sides["complex"] for g in str(cx).split("_") if g in gidx]
    gs.with_residuals(X, codes, np.concatenate([np.asarray(cand, np.int64)] + [np.asarray(g, np.int64) for g in groups.values()]
                                               + [np.asarray(g, np.int64) for g in groups_u.values()]), n_cells_corr, seed)
    res = side_tests(gs, sides, groups, matched=True)
    res_u = side_tests(gs, sides, groups_u, matched=False)
    res = res.assign(p_unmatched=res_u.p.to_numpy())
    key = pd.MultiIndex.from_frame(res[["complex", "type"]])

    def at(side, col):
        return res[col].to_numpy()[key.get_indexer(pd.MultiIndex.from_frame(side))]

    out = pd.DataFrame(index=rows.index)
    out["gate_b2"] = at(lig, "gate") & at(rec, "gate")
    out["p_L_b2"], out["p_R_b2"] = at(lig, "p"), at(rec, "p")
    out["p_b2"] = np.maximum(out.p_L_b2, out.p_R_b2)
    out["p_b2_uncorrected"] = np.maximum(at(lig, "p_uncorrected"), at(rec, "p_uncorrected"))
    out["p_b2_unmatched"] = np.maximum(at(lig, "p_unmatched"), at(rec, "p_unmatched"))
    zj = (at(lig, "z") + at(rec, "z")) / np.sqrt(2)
    out["p_b2_joint"] = np.where(out.gate_b2 & np.isfinite(zj), sstats.norm.sf(np.nan_to_num(zj)), 1.0)
    out["limiting_partner_b2"] = np.where(out.p_L_b2 >= out.p_R_b2, "ligand", "receptor")
    out.attrs["sides"] = res          # one row per (complex, type): the validity check reads side-level p
    return out
