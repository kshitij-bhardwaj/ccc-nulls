"""Amendment 1, A2: downstream-response truths. Does the receiving type t show the transcriptional
footprint of the call's signal (on genes other than L and R)?

Footprints:
- CytoSig (Jiang et al. 2021) ``signature.centroid``: 43 cytokine response signatures, human;
- PROGENy (Schubert et al. 2018) via OmniPath: 14 pathways, top 500 genes per pathway by p-value.
  An L-R is in a pathway, with a sign, when every L and R subunit has a weight of that sign there
  (LIANA+'s ``generate_lr_geneset`` rule).
Activity of a footprint in type t: the slope t-value of a univariate linear model (decoupler's
``ulm``) of t's contrast profile (mean log-expression in t minus the mean over the other types) on
the footprint weights, over the resource's genes, with the call's L and R genes removed. Per
footprint the activities are z-scored over the receiving types; t responds if z >= 1.645.
Mouse data is moved to human genes by HCOP one-to-one orthologs first.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import sparse

from ccc_nulls import datasets

Z_RESPONDS = 1.645
Z_SENSITIVITY = 1.0
PROGENY_TOP = 500

# CytoSig signature -> the ligand genes that carry it (NO is not a gene product: skipped)
CYTOKINE_GENES = {
    "Activin A": ["INHBA"], "BDNF": ["BDNF"], "BMP2": ["BMP2"], "BMP4": ["BMP4"], "BMP6": ["BMP6"],
    "CD40L": ["CD40LG"], "CXCL12": ["CXCL12"], "EGF": ["EGF"], "FGF2": ["FGF2"], "GCSF": ["CSF3"],
    "GDF11": ["GDF11"], "GMCSF": ["CSF2"], "HGF": ["HGF"],
    "IFN1": ["IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNA10", "IFNA13", "IFNA14",
             "IFNA16", "IFNA17", "IFNA21", "IFNB1", "IFNE", "IFNK", "IFNW1"],
    "IFNG": ["IFNG"], "IFNL": ["IFNL1", "IFNL2", "IFNL3", "IFNL4"], "IL10": ["IL10"], "IL12": ["IL12A", "IL12B"],
    "IL13": ["IL13"], "IL15": ["IL15"], "IL17A": ["IL17A"], "IL1A": ["IL1A"], "IL1B": ["IL1B"], "IL2": ["IL2"],
    "IL21": ["IL21"], "IL22": ["IL22"], "IL27": ["IL27", "EBI3"], "IL3": ["IL3"], "IL36": ["IL36A", "IL36B", "IL36G"],
    "IL4": ["IL4"], "IL6": ["IL6"], "LIF": ["LIF"], "LTA": ["LTA"], "MCSF": ["CSF1"], "OSM": ["OSM"],
    "TGFB1": ["TGFB1"], "TGFB3": ["TGFB3"], "TNFA": ["TNF"], "TRAIL": ["TNFSF10"], "TWEAK": ["TNFSF12"],
    "VEGFA": ["VEGFA"], "WNT3A": ["WNT3A"],
}


def data_dir():
    return datasets.data_dir() / "downstream"


@lru_cache(None)
def cytosig() -> pd.DataFrame:
    """genes x cytokines response signature (human symbols)."""
    return pd.read_csv(data_dir() / "cytosig_signature.centroid", sep="\t", index_col=0)


@lru_cache(None)
def progeny(top: int = PROGENY_TOP) -> pd.DataFrame:
    """genes x pathways weights (human), 0 outside a pathway's top genes by p-value."""
    p = pd.read_csv(data_dir() / "progeny_omnipath_human.tsv", sep="\t")
    w = p.pivot_table(index=["record_id", "genesymbol"], columns="label", values="value", aggfunc="first").reset_index()
    w["weight"], w["p_value"] = w.weight.astype(float), w.p_value.astype(float)
    w = w.sort_values(["pathway", "p_value"]).groupby("pathway").head(top)
    return w.pivot_table(index="genesymbol", columns="pathway", values="weight", aggfunc="first").fillna(0.0)


def human_rows(rows: pd.DataFrame, species: str) -> pd.DataFrame:
    """Rows with ligand/receptor subunits as human symbols (mouse CellChat rows carry mouse ones)."""
    if species == "human":
        return rows
    from .data import orthologs
    m = orthologs(species)
    mp = dict(zip(m.other_symbol, m.human_symbol))

    def conv(cx):
        out = [mp.get(g, g if g.isupper() else None) for g in str(cx).split("_")]
        return None if any(g is None for g in out) else "_".join(sorted(out))
    rows = rows.copy()
    for c in ("ligand_complex", "receptor_complex"):
        u = pd.Series(rows[c].unique())
        rows[c] = rows[c].map(dict(zip(u, u.map(conv))))
    return rows


def contrast_profiles(adata, groupby: str) -> pd.DataFrame:
    """genes x types: mean log-expression in the type minus the mean over the other types."""
    X = sparse.csr_matrix(adata.X)
    labels = adata.obs[groupby].astype(str).to_numpy()
    types = np.unique(labels)
    onehot = sparse.csr_matrix((np.ones(len(labels)), (np.searchsorted(types, labels), np.arange(len(labels)))),
                               shape=(len(types), len(labels)))
    mean = (onehot @ X).toarray() / np.asarray(onehot.sum(1))
    other = (mean.sum(0, keepdims=True) - mean) / (len(types) - 1)
    return pd.DataFrame((mean - other).T, index=adata.var_names.astype(str), columns=types)


def ulm(Y: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Slope t-value of y ~ a + b w for every column y of Y (genes x types)."""
    n = len(w)
    wc = w - w.mean()
    sxx = (wc ** 2).sum()
    if n < 3 or sxx == 0:
        return np.full(Y.shape[1], np.nan)
    Yc = Y - Y.mean(0)
    b = wc @ Yc / sxx
    resid = Yc - np.outer(wc, b)
    se = np.sqrt((resid ** 2).sum(0) / (n - 2) / sxx)
    with np.errstate(divide="ignore", invalid="ignore"):
        return b / se


class Footprints:
    """Activities of one resource's footprints in every type of one dataset, cached per (footprint, removed genes)."""

    def __init__(self, profiles: pd.DataFrame, weights: pd.DataFrame):
        genes = profiles.index.intersection(weights.index)
        self.genes = pd.Index(genes)
        self.Y = profiles.loc[genes].to_numpy(float)
        self.W = weights.loc[genes]
        self.types = profiles.columns.to_numpy()
        self._cache = {}

    def z(self, footprint: str, drop: tuple) -> pd.Series:
        key = (footprint, drop)
        if key not in self._cache:
            keep = ~self.genes.isin(drop)
            a = ulm(self.Y[keep], self.W[footprint].to_numpy(float)[keep])
            z = (a - np.nanmean(a)) / np.nanstd(a) if np.nanstd(a) > 0 else np.full_like(a, np.nan)
            self._cache[key] = (pd.Series(a, index=self.types), pd.Series(z, index=self.types))
        return self._cache[key]


def lr_pathways(lr: pd.DataFrame, weights: pd.DataFrame) -> pd.DataFrame:
    """(ligand_complex, receptor_complex, pathway, sign): every subunit weighted with one sign there."""
    out = []
    for L, R in lr.drop_duplicates().itertuples(index=False):
        genes = str(L).split("_") + str(R).split("_")
        if not all(g in weights.index for g in genes):
            continue
        W = weights.loc[genes]
        for pw in W.columns:
            s = np.sign(W[pw].to_numpy())
            if (s != 0).all() and (s == s[0]).all():
                out.append((L, R, pw, int(s[0])))
    return pd.DataFrame(out, columns=["ligand_complex", "receptor_complex", "footprint", "sign"])


def truth_rows(rows: pd.DataFrame, profiles: pd.DataFrame, resource: str) -> pd.DataFrame:
    """Null A's calls with s != t that the resource covers, one line per (call, footprint), with the
    receiver's response (z >= 1.645; z >= 1.0 as sensitivity)."""
    r = rows[(rows.source != rows.target) & rows.call_null_a.to_numpy(bool)
             & rows.ligand_complex.notna() & rows.receptor_complex.notna()].copy()
    if resource == "cytosig":
        weights = cytosig()
        gene2cyt = {g: c for c, gs in CYTOKINE_GENES.items() for g in gs if c in weights.columns}
        cyt = r.ligand_complex.map(lambda cx: next((gene2cyt[g] for g in str(cx).split("_") if g in gene2cyt), None))
        r = r.assign(footprint=cyt, sign=1)[cyt.notna()]
    else:
        weights = progeny()
        m = lr_pathways(r[["ligand_complex", "receptor_complex"]], weights)
        r = r.merge(m, on=["ligand_complex", "receptor_complex"])
    fp = Footprints(profiles, weights)
    act, z = [], []
    for row in r[["footprint", "ligand_complex", "receptor_complex", "target", "sign"]].itertuples(index=False):
        drop = tuple(sorted(set(str(row.ligand_complex).split("_")) | set(str(row.receptor_complex).split("_"))))
        a, zz = fp.z(row.footprint, drop)
        a, zz = a * row.sign, zz * row.sign
        act.append(a.get(row.target, np.nan))
        z.append(zz.get(row.target, np.nan))
    r["activity"], r["z"] = act, z
    r = r[np.isfinite(r.z)]
    r["responds"] = r.z >= Z_RESPONDS
    r["responds_z1"] = r.z >= Z_SENSITIVITY
    r["unit"] = r.footprint.astype(str) + " -> " + r.target.astype(str)
    return r.reset_index(drop=True)


def odds_ratio(r: pd.DataFrame, kept: np.ndarray, dropped: np.ndarray, col: str = "responds",
               n_boot: int = 2000, seed: int = 0) -> dict:
    """Odds of response, kept vs dropped calls; bootstrap over truth units (footprint, receiver)."""
    d = r.assign(k=kept.astype(int), dr=dropped.astype(int), y=r[col].astype(int))
    d = d[(d.k + d.dr) > 0]
    g = d.groupby("unit").apply(lambda x: pd.Series(dict(a=(x.k * x.y).sum(), b=(x.k * (1 - x.y)).sum(),
                                                          c=(x.dr * x.y).sum(), e=(x.dr * (1 - x.y)).sum())),
                                include_groups=False)
    K = g.to_numpy(float)

    def orr(k):
        a, b, c, e = k + (0.5 if (k == 0).any() else 0)
        return (a / b) / (c / e)
    tot = K.sum(0)
    if len(K) == 0 or tot[:2].sum() == 0 or tot[2:].sum() == 0:
        return dict(odds_ratio=None, ci95=None, kept_calls=int(tot[:2].sum()) if len(K) else 0,
                    dropped_calls=int(tot[2:].sum()) if len(K) else 0, n_units=int(len(K)))
    rng = np.random.default_rng(seed)
    boot = np.array([orr(K[rng.integers(0, len(K), len(K))].sum(0)) for _ in range(n_boot)])
    return dict(odds_ratio=float(orr(tot)), ci95=[float(x) for x in np.percentile(boot, [2.5, 97.5])],
                kept_calls=int(tot[0] + tot[1]), kept_responding=int(tot[0]),
                dropped_calls=int(tot[2] + tot[3]), dropped_responding=int(tot[2]),
                frac_responding_kept=float(tot[0] / max(tot[0] + tot[1], 1)),
                frac_responding_dropped=float(tot[2] / max(tot[2] + tot[3], 1)), n_units=int(len(K)))


def endpoints(r: pd.DataFrame) -> dict:
    """D1 (kept vs dropped) and D2 (kept vs receptor-limited / ligand-limited dropped), z >= 1.645,
    plus D1 at the sensitivity threshold."""
    kept = r.call_null_b.to_numpy(bool)
    dropped = ~kept
    rec = dropped & (r.limiting_partner.to_numpy() == "receptor")
    lig = dropped & (r.limiting_partner.to_numpy() == "ligand")
    return dict(n_lines=int(len(r)), n_units=int(r.unit.nunique()), n_footprints=int(r.footprint.nunique()),
                responding=int(r.responds.sum()),
                D1=odds_ratio(r, kept, dropped), D1_z1=odds_ratio(r, kept, dropped, col="responds_z1"),
                D2_receptor=odds_ratio(r, kept, rec), D2_ligand=odds_ratio(r, kept, lig))
