"""A single-cell spatial truth for TNBC from 10x Xenium human breast cancer (Janesick et al., Nat Commun 2023).

NOT pre-registered: an addition. Visium spots mix several cells, so TNBC's existing spatial truth is inferred
by deconvolution; Xenium resolves single cells. It is a different patient and not triple negative, so TNBC's
29 labels and Xenium's 20 are both mapped to 10 broad classes present in both (``TNBC_CLASS``, ``XENIUM_CLASS``);
labels without a counterpart are left out.

Data (data/xenium_breast/):
- cells_rep{1,2}.csv.gz: Xenium_FFPE_Human_Breast_Cancer_Rep{1,2}_cells.csv.gz (cf.10xgenomics.com, samples/xenium/1.0.1)
- GSE243275_Barcode_Cell_Type_Matrices.xlsx (GEO GSE243275): sheets "Xenium R{1,2} Fig1-5 (supervised)"

Truth: the hypothalamus MERFISH method (``hypothalamus.neighbourhood_truth``), with the two tissue sections as
slices: 10-nearest-neighbour graph on centroids, edges between classes s and t counted (both directions), z against
200 label permutations within the section; mean z over sections where both classes have >= 10 cells;
colocalised if >= 1.645. Unlabelled and hybrid cells stay in the graph but are not counted. The class-level truth
is expanded to TNBC label pairs (s, t) whose classes differ.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from ccc_nulls import datasets

K_NN, N_PERM, MIN_CELLS, Z_COLOCALISED = 10, 200, 10, 1.645

XENIUM_CLASS = {
    "Invasive_Tumor": "Tumour", "Prolif_Invasive_Tumor": "Tumour", "DCIS_1": "Tumour", "DCIS_2": "Tumour",
    "Myoepi_ACTA2+": "Myoepithelial", "Myoepi_KRT15+": "Myoepithelial", "Stromal": "Fibroblast",
    "Perivascular-Like": "Perivascular", "Endothelial": "Endothelial", "Macrophages_1": "Myeloid",
    "Macrophages_2": "Myeloid", "IRF7+_DCs": "Dendritic", "LAMP3+_DCs": "Dendritic", "CD4+_T_Cells": "CD4 T",
    "CD8+_T_Cells": "CD8 T", "B_Cells": "B",
}
TNBC_CLASS = {
    "Cancer.Basal.SC": "Tumour", "Cancer.Her2.SC": "Tumour", "Cancer.LumA.SC": "Tumour", "Cancer.LumB.SC": "Tumour",
    "Cancer.Cycling": "Tumour", "Myoepithelial": "Myoepithelial", "CAFs.MSC.iCAF.like": "Fibroblast",
    "CAFs.myCAF.like": "Fibroblast", "PVL.Differentiated": "Perivascular", "PVL.Immature": "Perivascular",
    "Cycling.PVL": "Perivascular", "Endothelial.ACKR1": "Endothelial", "Endothelial.CXCL12": "Endothelial",
    "Endothelial.Lymphatic.LYVE1": "Endothelial", "Endothelial.RGS5": "Endothelial", "Macrophage": "Myeloid",
    "Monocyte": "Myeloid", "Cycling.Myeloid": "Myeloid", "DCs": "Dendritic", "T.cells.CD4": "CD4 T",
    "T.cells.CD8": "CD8 T", "B.cells.Memory": "B", "B.cells.Naive": "B",
}


def root():
    return datasets.data_dir() / "xenium_breast"


def sections() -> dict:
    """section -> DataFrame(x, y, cls) with cls NaN for cells not counted."""
    x = pd.ExcelFile(root() / "GSE243275_Barcode_Cell_Type_Matrices.xlsx")
    out = {}
    for rep in (1, 2):
        lab = x.parse(f"Xenium R{rep} Fig1-5 (supervised)").set_axis(["cell_id", "label"], axis=1)
        cells = pd.read_csv(root() / f"cells_rep{rep}.csv.gz", usecols=["cell_id", "x_centroid", "y_centroid"])
        d = cells.merge(lab, on="cell_id", how="left")
        d["cls"] = d.label.map(XENIUM_CLASS)
        out[f"rep{rep}"] = d
    return out


def class_truth(log=print) -> tuple[pd.DataFrame, dict]:
    from sklearn.neighbors import NearestNeighbors

    secs = sections()
    classes = np.array(sorted(set(XENIUM_CLASS.values())))
    K = len(classes)
    rng = np.random.default_rng(0)
    zs, info = [], {"sections": {}}
    for name, d in secs.items():
        xy = d[["x_centroid", "y_centroid"]].to_numpy(float)
        lab = pd.Categorical(d.cls, categories=classes).codes
        nn = NearestNeighbors(n_neighbors=K_NN + 1).fit(xy).kneighbors(xy, return_distance=False)[:, 1:]
        i, j = np.repeat(np.arange(len(xy)), K_NN), nn.ravel()

        def counts(l):
            a, b = l[i], l[j]
            ok = (a >= 0) & (b >= 0)
            c = np.bincount(a[ok] * K + b[ok], minlength=K * K).reshape(K, K).astype(float)
            return c + c.T
        obs = counts(lab)
        perm = np.stack([counts(rng.permutation(lab)) for _ in range(N_PERM)])
        sd = perm.std(0)
        z = np.divide(obs - perm.mean(0), sd, out=np.full((K, K), np.nan), where=sd > 0)
        n = np.bincount(lab[lab >= 0], minlength=K)
        z[(n < MIN_CELLS)[:, None] | (n < MIN_CELLS)[None, :]] = np.nan
        zs.append(z)
        info["sections"][name] = dict(cells=int(len(d)), labelled=int((d.label.notna()).sum()),
                                      counted=int((lab >= 0).sum()), cells_per_class=dict(zip(classes, n.tolist())))
    Z = np.stack(zs)
    zmean, nsec = np.nanmean(Z, 0), np.isfinite(Z).sum(0)
    s, t = np.meshgrid(np.arange(K), np.arange(K), indexing="ij")
    truth = pd.DataFrame({"source_class": classes[s.ravel()], "target_class": classes[t.ravel()],
                          "z_mean": zmean.ravel(), "n_sections": nsec.ravel()})
    truth = truth[(truth.source_class != truth.target_class) & (truth.n_sections > 0)].reset_index(drop=True)
    truth["colocalised"] = truth.z_mean >= Z_COLOCALISED
    info.update(n_pairs=int(len(truth)), n_colocalised=int(truth.colocalised.sum()), k_nn=K_NN, n_perm=N_PERM,
                z_threshold=Z_COLOCALISED, classes=list(classes))
    log(f"  Xenium truth: {info['n_colocalised']} of {info['n_pairs']} ordered class pairs colocalised")
    return truth, info


def tnbc_truth(out_dir, log=print) -> pd.DataFrame:
    """The class truth expanded to TNBC label pairs (labels as in the benchmark's data: '|' -> '/')."""
    cls, info = class_truth(log)
    lab = pd.Series(TNBC_CLASS)
    pairs = pd.DataFrame([(a, b) for a in lab.index for b in lab.index if lab[a] != lab[b]], columns=["source", "target"])
    pairs["source_class"], pairs["target_class"] = lab[pairs.source].to_numpy(), lab[pairs.target].to_numpy()
    truth = pairs.merge(cls, on=["source_class", "target_class"])
    cls.to_csv(out_dir / "truth_xenium_breast_classes.csv", index=False)
    truth.to_csv(out_dir / "truth_tnbc_wu2021_xenium.csv", index=False)
    (out_dir / "truth_tnbc_wu2021_xenium.json").write_text(json.dumps(info, indent=1))
    return truth
