"""Spatial ground truth for (sender, receiver) cell-type pairs: which pairs are colocalised.

TNBC (Wu et al. 2021): the study's own 10x Visium slides (Zenodo 4739739), the 4 TNBC ones, each
spot deconvolved onto the scRNA atlas's cell types, then the procedure of Dimitrov et al. 2022
(Nat Commun 13:3224), where LIANA was benchmarked against these slides:
- deconvolution: non-negative least squares of each spot's CP10k on the cell types' mean CP10k
  profiles over marker genes (the regression step of SPOTlight, which that paper used);
- per slide, the Pearson correlation of cell-type proportions across spots, as z-scores over all
  s != t pairs; colocalised = mean z over slides >= 1.645.
Allen brain: Open Problems' own source-target truth (adata.uns["ccc_target"]), used as given.

    python -m benchmarks.liana_nulls.spatial      # writes results/liana_nulls/truth_tnbc_spatial.csv
"""
from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio
from scipy import sparse
from scipy.optimize import nnls

from ccc_nulls import datasets

TNBC_SLIDES = ["1142243F", "1160920F", "CID4465", "CID44971"]   # subtype TNBC in the slides' metadata
Z_COLOCALISED = 1.645
N_MARKERS = 50


def visium_dir() -> Path:
    return Path(datasets.data_dir()) / "wu2021_visium"


def _open(path: Path):
    """The Zenodo files are named .gz but are plain text; open either kind."""
    with open(path, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rb") if gz else open(path, "rb")


def read_slide(slide: str):
    """Spots x genes counts (CSR), gene symbols, and the spot metadata (pathology, subtype)."""
    d = visium_dir() / "filtered_count_matrices" / f"{slide}_filtered_count_matrix"
    m = sio.mmread(_open(d / "matrix.mtx.gz")).T.tocsr()
    genes = pd.read_csv(_open(d / "features.tsv.gz"), sep="\t", header=None)
    barcodes = pd.read_csv(_open(d / "barcodes.tsv.gz"), sep="\t", header=None)[0].to_numpy()
    symbols = genes[1 if genes.shape[1] > 1 else 0].astype(str).to_numpy()
    meta = pd.read_csv(visium_dir() / "metadata" / f"{slide}_metadata.csv", index_col=0)
    return m, symbols, barcodes, meta.reindex(barcodes)


def reference_profiles(adata, groupby: str, n_markers: int = N_MARKERS):
    """Cell-type mean CP10k over the union of each type's top markers (mean log-expression in the type
    minus the mean over the other types, expressed in >= 20% of the type's cells)."""
    X = sparse.csr_matrix(adata.X)
    labels = adata.obs[groupby].astype(str).to_numpy()
    types = np.unique(labels)
    onehot = sparse.csr_matrix((np.ones(len(labels)), (np.searchsorted(types, labels), np.arange(len(labels)))),
                               shape=(len(types), len(labels)))
    n = np.asarray(onehot.sum(1)).ravel()[:, None]
    mean_log = (onehot @ X).toarray() / n
    frac = (onehot @ (X > 0).astype(np.float32)).toarray() / n
    mean_cp10k = (onehot @ X.expm1()).toarray() / n
    markers = set()
    for k in range(len(types)):
        other = np.delete(mean_log, k, 0).mean(0)
        score = np.where(frac[k] >= 0.2, mean_log[k] - other, -np.inf)
        markers.update(np.argsort(score)[::-1][:n_markers][np.isfinite(np.sort(score)[::-1][:n_markers])].tolist())
    markers = np.array(sorted(markers))
    return pd.DataFrame(mean_cp10k[:, markers].T, index=adata.var_names[markers].astype(str), columns=types)


def deconvolve(counts, symbols, profiles: pd.DataFrame) -> np.ndarray:
    """Spots x types proportions: NNLS of each spot's CP10k (marker genes) on the type profiles."""
    sym = pd.Index(symbols)
    first = ~sym.duplicated()
    pos = pd.Series(np.arange(len(sym))[first], index=sym[first])
    genes = profiles.index[profiles.index.isin(pos.index)]
    lib = np.asarray(counts.sum(1)).ravel()
    Y = np.asarray(counts[:, pos[genes].to_numpy()].todense()) / np.maximum(lib, 1)[:, None] * 1e4
    M = profiles.loc[genes].to_numpy()
    W = np.array([nnls(M, y)[0] for y in Y])
    s = W.sum(1, keepdims=True)
    return np.divide(W, s, out=np.zeros_like(W), where=s > 0), len(genes)


def colocalisation(props: np.ndarray, types) -> pd.DataFrame:
    """Per slide: Pearson correlation of proportions across spots, z-scored over the s != t pairs."""
    C = pd.DataFrame(props, columns=types).corr().to_numpy()
    off = ~np.eye(len(types), dtype=bool) & np.isfinite(C)
    z = np.full(C.shape, np.nan)
    z[off] = (C[off] - C[off].mean()) / C[off].std()
    s, t = np.meshgrid(np.arange(len(types)), np.arange(len(types)), indexing="ij")
    return pd.DataFrame({"source": np.asarray(types)[s.ravel()], "target": np.asarray(types)[t.ravel()],
                         "r": C.ravel(), "z": z.ravel()})


def tnbc_truth(adata, groupby: str, slides=TNBC_SLIDES, log=print) -> tuple[pd.DataFrame, dict]:
    profiles = reference_profiles(adata, groupby)
    types = profiles.columns.to_numpy()
    per, info = [], {"n_marker_genes": int(len(profiles)), "slides": {}}
    for sl in slides:
        counts, symbols, _, meta = read_slide(sl)
        props, n_used = deconvolve(counts, symbols, profiles)
        per.append(colocalisation(props, types).assign(slide=sl))
        top = pd.Series(props.mean(0), index=types).nlargest(5).round(3).to_dict()
        info["slides"][sl] = dict(spots=int(counts.shape[0]), marker_genes_used=int(n_used),
                                  subtype=str(meta.subtype.dropna().iloc[0]) if "subtype" in meta else None,
                                  mean_proportion_top5=top)
        log(f"    {sl}: {counts.shape[0]:,} spots, {n_used} markers; top types {top}")
    z = pd.concat(per)
    truth = z.groupby(["source", "target"], sort=False).agg(r_mean=("r", "mean"), z_mean=("z", "mean"),
                                                             n_slides=("z", "count")).reset_index()
    truth = truth[truth.source != truth.target]
    truth["colocalised"] = truth.z_mean >= Z_COLOCALISED
    info.update(n_pairs=int(len(truth)), n_colocalised=int(truth.colocalised.sum()), z_threshold=Z_COLOCALISED)
    return truth.reset_index(drop=True), info


def brain_truth(adata) -> pd.DataFrame:
    t = adata.uns["ccc_target"].rename(columns={"response": "colocalised"}).copy()
    t["colocalised"] = t.colocalised.astype(int).astype(bool)
    return t[t.source != t.target].reset_index(drop=True)


if __name__ == "__main__":
    import json

    from .data import load
    adata, groupby, _ = load("tnbc_wu2021")
    truth, info = tnbc_truth(adata, groupby)
    out = Path(__file__).resolve().parents[1] / "results" / "liana_nulls"
    out.mkdir(parents=True, exist_ok=True)
    truth.to_csv(out / "truth_tnbc_spatial.csv", index=False)
    (out / "truth_tnbc_spatial.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))
