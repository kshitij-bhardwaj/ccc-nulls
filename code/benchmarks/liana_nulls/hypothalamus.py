"""Amendment 1, A1: the mouse hypothalamic preoptic region of Moffitt et al., Science 2018.

scRNA-seq (GEO GSE113576, 10x) carries the nulls; MERFISH (animal 1, 12 slices, 161 genes, the
study's own ``Cell_class``) is the spatial truth only.

scRNA cell classes come from the MERFISH classes by transfer (the study's per-cell scRNA labels are
not downloadable here): log CP10k -> 2,000 HVGs -> 50 PCs -> k-means (k = 40, seed 0); each
cluster's mean over the shared panel genes and each MERFISH class's mean are z-scored per gene; a
cluster takes the class of highest Pearson r if r >= 0.3, else its cells are left out.

Truth: per slice, 10-nearest-neighbour graph on the centroids, edges between classes s and t
counted (both directions), z against 200 label permutations within the slice; mean z over the slices
where both classes have >= 10 cells; colocalised if >= 1.645.
"""
from __future__ import annotations

import gzip
import json

import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse

from ccc_nulls import datasets

NAME = "moffitt2018_hypothalamus"
MERGE = {"OD Mature 1": "OD Mature", "OD Mature 2": "OD Mature", "OD Mature 3": "OD Mature", "OD Mature 4": "OD Mature",
         "OD Immature 1": "OD Immature", "OD Immature 2": "OD Immature",
         "Endothelial 1": "Endothelial", "Endothelial 2": "Endothelial", "Endothelial 3": "Endothelial"}
MARKERS = ["Gad1", "Slc17a6", "Aqp4", "Mlc1", "Mbp", "Ermn", "Pdgfra", "Slco1a4", "Cd24a", "Selplg", "Rgs5"]
K_CLUSTERS, R_MIN = 40, 0.3
K_NN, N_PERM, MIN_CELLS, Z_COLOCALISED = 10, 200, 10, 1.645


def root():
    return datasets.data_dir() / "moffitt2018"


def merfish():
    import anndata as ad
    a = ad.read_h5ad(root() / "merfish" / "merfish.h5ad")
    cls = a.obs.Cell_class.astype(str).replace(MERGE)
    a.obs["cell_class"] = cls.where(cls != "Ambiguous")
    return a


def read_scrna():
    """cells x genes counts (CSR) and gene symbols from the GEO 10x files."""
    d = root() / "scrna"
    m = sio.mmread(gzip.open(d / "GSE113576_matrix.mtx.gz")).T.tocsr()
    genes = pd.read_csv(d / "GSE113576_genes.tsv.gz", sep="\t", header=None)
    bars = pd.read_csv(d / "GSE113576_barcodes.tsv.gz", sep="\t", header=None)[0].astype(str).to_numpy()
    return m, genes[genes.shape[1] - 1].astype(str).to_numpy(), bars


def _zrows(M):
    sd = M.std(0, keepdims=True)
    return np.divide(M - M.mean(0, keepdims=True), sd, out=np.zeros_like(M), where=sd > 0)


def build(log=print):
    """The scRNA AnnData (X = ln(1 + CP10k), obs cell_class) with transferred labels; cached in data/."""
    import anndata as ad
    import scanpy as sc
    from sklearn.cluster import KMeans

    f = root() / f"{NAME}.h5ad"
    if f.exists():
        return ad.read_h5ad(f)
    counts, symbols, bars = read_scrna()
    first = ~pd.Index(symbols).duplicated()
    a = ad.AnnData(sparse.csr_matrix(counts[:, first], dtype=np.float32), obs=pd.DataFrame(index=bars),
                   var=pd.DataFrame(index=symbols[first]))
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    h = a.copy()
    sc.pp.highly_variable_genes(h, n_top_genes=2000)
    h = h[:, h.var.highly_variable].copy()
    sc.pp.scale(h, max_value=10)
    sc.tl.pca(h, n_comps=50, random_state=0)
    cl = KMeans(K_CLUSTERS, random_state=0, n_init=10).fit_predict(h.obsm["X_pca"])
    del h

    m = merfish()
    shared = [g for g in m.var_names if g in a.var_names and not g.startswith("Blank")]
    mm = m[m.obs.cell_class.notna(), shared]
    Xm = np.log1p(mm.X.toarray() if sparse.issparse(mm.X) else np.asarray(mm.X))
    classes = np.array(sorted(mm.obs.cell_class.unique()))
    C = np.vstack([Xm[mm.obs.cell_class.to_numpy() == c].mean(0) for c in classes])          # classes x genes
    Xs = a[:, shared].X
    S = np.vstack([np.asarray(Xs[cl == k].mean(0)).ravel() for k in range(K_CLUSTERS)])     # clusters x genes
    Cz, Sz = _zrows(C), _zrows(S)
    R = np.corrcoef(np.vstack([Sz, Cz]))[:K_CLUSTERS, K_CLUSTERS:]                          # clusters x classes
    best, r = R.argmax(1), R.max(1)
    lab = np.where(r >= R_MIN, classes[best], None)
    a.obs["cluster"] = cl
    a.obs["cell_class"] = pd.Categorical(lab[cl])
    a.obs["transfer_r"] = r[cl]
    present = [c for c in classes if (a.obs.cell_class == c).any()]       # a class may receive no cluster
    marker = {g: {c: float(np.asarray(a[a.obs.cell_class == c, g].X.mean())) for c in present}
              for g in MARKERS if g in a.var_names}
    info = dict(n_cells=int(a.n_obs), n_genes=int(a.n_vars), shared_panel_genes=len(shared),
                clusters={int(k): dict(cell_class=None if lab[k] is None else str(lab[k]), r=round(float(r[k]), 3),
                                       n=int((cl == k).sum())) for k in range(K_CLUSTERS)},
                cells_per_class=a.obs.cell_class.value_counts().to_dict(),
                classes_without_cells=[str(c) for c in classes if c not in present], unassigned=int(a.obs.cell_class.isna().sum()),
                marker_mean_log_cp10k=marker)
    a = a[a.obs.cell_class.notna()].copy()
    a.obs["cell_class"] = a.obs.cell_class.astype(str).astype("category")
    a.uns["transfer"] = json.dumps(info)
    a.write_h5ad(f, compression="gzip")
    log(f"  hypothalamus scRNA: {a.n_obs:,} cells labelled, {info['unassigned']:,} unassigned; {info['cells_per_class']}")
    return a


def neighbourhood_truth(log=print) -> tuple[pd.DataFrame, dict]:
    from sklearn.neighbors import NearestNeighbors

    m = merfish()
    obs = m.obs
    classes = np.array(sorted(obs.cell_class.dropna().unique()))
    K = len(classes)
    rng = np.random.default_rng(0)
    zs, info = [], {"slices": {}}
    for sl, g in obs.groupby("Bregma", observed=True):
        xy = g[["Centroid_X", "Centroid_Y"]].to_numpy(float)
        lab = pd.Categorical(g.cell_class, categories=classes).codes      # -1 = Ambiguous: in the graph, not counted
        nn = NearestNeighbors(n_neighbors=K_NN + 1).fit(xy).kneighbors(xy, return_distance=False)[:, 1:]
        i, j = np.repeat(np.arange(len(xy)), K_NN), nn.ravel()

        def counts(l):
            a, b = l[i], l[j]
            ok = (a >= 0) & (b >= 0)
            c = np.bincount(a[ok] * K + b[ok], minlength=K * K).reshape(K, K).astype(float)
            return c + c.T
        obs_c = counts(lab)
        perm = np.stack([counts(rng.permutation(lab)) for _ in range(N_PERM)])
        sd = perm.std(0)
        z = np.divide(obs_c - perm.mean(0), sd, out=np.full((K, K), np.nan), where=sd > 0)
        n = np.bincount(lab[lab >= 0], minlength=K)
        z[(n < MIN_CELLS)[:, None] | (n < MIN_CELLS)[None, :]] = np.nan
        zs.append(z)
        info["slices"][str(round(float(sl), 1))] = dict(cells=int(len(g)), cells_per_class=dict(zip(classes, n.tolist())))
    Z = np.stack(zs)
    zmean, nsl = np.nanmean(Z, 0), np.isfinite(Z).sum(0)
    s, t = np.meshgrid(np.arange(K), np.arange(K), indexing="ij")
    truth = pd.DataFrame({"source": classes[s.ravel()], "target": classes[t.ravel()], "z_mean": zmean.ravel(),
                          "n_slices": nsl.ravel()})
    truth = truth[(truth.source != truth.target) & (truth.n_slices > 0)].reset_index(drop=True)
    truth["colocalised"] = truth.z_mean >= Z_COLOCALISED
    info.update(n_pairs=int(len(truth)), n_colocalised=int(truth.colocalised.sum()), k_nn=K_NN, n_perm=N_PERM,
                z_threshold=Z_COLOCALISED)
    log(f"  hypothalamus truth: {info['n_colocalised']} of {info['n_pairs']} ordered class pairs colocalised")
    return truth, info
