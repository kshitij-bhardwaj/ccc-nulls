"""Dataset registry, downloads and loaders.

Files live in ``<repo>/data/`` by default; set ``CCC_NULLS_DATA`` to use another folder.

    python -m ccc_nulls.datasets tnbc_wu2021          # download + md5 check
    python -m ccc_nulls.datasets tabula_sapiens_mammary tabula_sapiens_liver
    >>> from ccc_nulls.datasets import load_openproblems, load_tabula_sapiens
    >>> adata, truth = load_openproblems("tnbc_wu2021")
    >>> adata = load_tabula_sapiens("tabula_sapiens_mammary")
"""
from __future__ import annotations

import hashlib
import os
import sys
import tempfile
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

S3 = "https://openproblems-data.s3.amazonaws.com/resources/datasets/openproblems_v1"
TABULA_SAPIENS = "https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5"

DATASETS = {
    "tnbc_wu2021": dict(
        role="evaluation: cytokine ground truth (ligand -> target cell type)",
        url=f"{S3}/tnbc_wu2021/log_cp10k/dataset.h5ad",
        filename="tnbc_wu2021_log_cp10k.h5ad",
        md5="0c1c558c4137571eb5e9351c8408dc04",
        size=1_269_138_232,
        organism="human",
        groupby="label",
        sample_key="orig.ident",
        truth_keys=["ligand", "target"],
        summary="Triple-negative breast cancer atlas (Wu et al. 2021): 42,512 cells, 29 cell types, "
                "10 patients. Open Problems ligand-target task; used in LIANA's benchmark.",
    ),
    "allen_brain_atlas": dict(
        role="evaluation: spatial ground truth (source -> target adjacency)",
        url=f"{S3}/allen_brain_atlas/log_cp10k/dataset.h5ad",
        filename="allen_brain_atlas_log_cp10k.h5ad",
        md5="89a5f0919ed3f2b3e03e7550ef5a20ef",
        size=2_049_800_441,
        organism="mouse",
        groupby="label",
        sample_key=None,
        truth_keys=["source", "target"],
        summary="Adult mouse primary visual cortex (GSE71585): 14,249 cells. Truth = adjacent cell types "
                "from deconvolved 10x Visium. Open Problems source-target task.",
    ),
    "tabula_sapiens": dict(
        role="calibration: no ground truth (permutation / fake-resource tests)",
        url=TABULA_SAPIENS,
        filename=None,
        md5=None,
        size=None,
        organism="human",
        groupby="cell_type",
        sample_key="donor_id",
        truth_keys=None,
        summary="Tabula Sapiens (CELLxGENE), one h5ad per tissue. The tissues in use are registered as "
                "tabula_sapiens_<tissue>; load them with load_tabula_sapiens().",
    ),
    # Tabula Sapiens tissues used for calibration (CELLxGENE files of the collection above)
    "tabula_sapiens_mammary": dict(
        role="calibration: no ground truth (full-permutation and known-null tests)",
        url="https://datasets.cellxgene.cziscience.com/33b25f68-a9d6-4643-9283-653e937ae608.h5ad",
        filename="tabula_sapiens_mammary.h5ad",
        md5="6261c8dd5785c52fb8a18d0974602597",
        size=1_334_205_970,
        organism="human",
        groupby="label",
        sample_key="donor_id",
        truth_keys=None,
        collection=TABULA_SAPIENS,
        title="Tabula Sapiens - Mammary",
        summary="Normal breast, the tissue closest to the TNBC atlas: 30,936 cells, 2 donors. Raw counts in "
                "raw/X; load_tabula_sapiens() normalises them to ln(1 + CP10k).",
    ),
    "tabula_sapiens_liver": dict(
        role="calibration: no ground truth (full-permutation and known-null tests)",
        url="https://datasets.cellxgene.cziscience.com/a4b3a49e-062b-4e3f-8915-02f40607651f.h5ad",
        filename="tabula_sapiens_liver.h5ad",
        md5="ccb220e15f2070c3c22c226fe89d5878",
        size=834_485_899,
        organism="human",
        groupby="label",
        sample_key="donor_id",
        truth_keys=None,
        collection=TABULA_SAPIENS,
        title="Tabula Sapiens - Liver",
        summary="A second organ of similar size: 22,214 cells, 4 donors. Raw counts in raw/X; "
                "load_tabula_sapiens() normalises them to ln(1 + CP10k).",
    ),
    "tabula_sapiens_lung": dict(
        role="held-out confirmation (benchmarks/liana_final/AMENDMENT_2.md): no ground truth",
        url="https://datasets.cellxgene.cziscience.com/fddf43a3-b7c8-49e9-ac0a-5e24b390849f.h5ad",
        filename="tabula_sapiens_lung.h5ad",
        md5="2efb556c5950390d8d42ca7828c3f3f3",
        size=3_196_847_981,
        organism="human",
        groupby="label",
        sample_key="donor_id",
        truth_keys=None,
        collection=TABULA_SAPIENS,
        title="Tabula Sapiens - Lung",
        summary="A third tissue, untouched until AMENDMENT_2's rule was frozen: 65,847 cells, 4 donors. Raw "
                "counts in raw/X; load_tabula_sapiens() normalises them to ln(1 + CP10k).",
    ),
}


def data_dir() -> Path:
    env = os.environ.get("CCC_NULLS_DATA")
    return Path(env) if env else Path(__file__).resolve().parents[2] / "data"


def path(name: str) -> Path:
    meta = DATASETS[name]
    if not meta["filename"]:
        raise ValueError(f"{name} has no single canonical file; see DATASETS['{name}']['summary']")
    return data_dir() / meta["filename"]


def md5sum(p: Path, chunk: int = 1 << 24) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _fetch_range(url: str, start: int, end: int, dest: Path, retries: int = 5) -> None:
    want = end - start + 1
    for _ in range(retries):
        req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
        try:
            with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
                while block := r.read(1 << 20):
                    f.write(block)
            if dest.stat().st_size == want:
                return
        except OSError:
            pass
    raise RuntimeError(f"failed to fetch bytes {start}-{end} of {url}")


def download(name: str, parts: int = 16, force: bool = False) -> Path:
    """Download with parallel byte ranges (a single stream from S3 can be very slow)."""
    meta = DATASETS[name]
    out = path(name)
    if out.exists() and not force:
        if meta["md5"] and md5sum(out) != meta["md5"]:
            raise RuntimeError(f"{out} exists but its md5 does not match; re-run with force=True")
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    size, step = meta["size"], -(-meta["size"] // parts)
    with tempfile.TemporaryDirectory(dir=out.parent) as tmp:
        chunks = [(i * step, min((i + 1) * step, size) - 1, Path(tmp) / f"part{i:02d}") for i in range(parts)]
        with ThreadPoolExecutor(parts) as ex:
            list(ex.map(lambda c: _fetch_range(meta["url"], *c), chunks))
        with open(out, "wb") as f:
            for _, _, p in chunks:
                f.write(p.read_bytes())
    if meta["md5"] and md5sum(out) != meta["md5"]:
        out.unlink()
        raise RuntimeError(f"md5 mismatch for {name}; file removed")
    return out


def load_h5ad(p, layer: str | None = "normalized"):
    """Read an h5ad and return a lean AnnData whose X is ``layer`` (log-normalised, CSR float32)."""
    import anndata as ad
    import numpy as np

    raw = ad.read_h5ad(p)
    X = raw.layers[layer] if layer else raw.X
    adata = ad.AnnData(X=X.astype(np.float32).tocsr(), obs=raw.obs.copy(), var=raw.var[[]].copy(),
                       uns={k: raw.uns[k] for k in ("ccc_target",) if k in raw.uns})
    return adata


def load_openproblems(name: str, layer: str = "normalized", min_cells: int = 5):
    """Load an Open Problems CCC dataset: (AnnData with X = log-CP10k, ground-truth DataFrame).

    Cell types with fewer than ``min_cells`` cells are dropped, as LIANA does.
    """
    meta = DATASETS[name]
    adata = load_h5ad(path(name), layer)
    vc = adata.obs[meta["groupby"]].value_counts()
    adata = adata[adata.obs[meta["groupby"]].isin(vc[vc >= min_cells].index)].copy()
    adata.obs[meta["groupby"]] = adata.obs[meta["groupby"]].astype(str).astype("category")
    truth = adata.uns["ccc_target"][meta["truth_keys"] + ["response"]].copy()
    truth["response"] = truth["response"].astype(int)
    return adata, truth.reset_index(drop=True)


def _h5_strings(f, key):
    """A string column of an h5ad dataframe, plain or categorical."""
    import h5py
    import numpy as np

    o = f[key]
    if isinstance(o, h5py.Group):                          # categorical: codes + categories
        return np.asarray(o["categories"].asstr()[:], dtype=object)[o["codes"][:]]
    return np.asarray(o.asstr()[:], dtype=object)


def load_tabula_sapiens(name: str, genes=None, min_cells: int = 5, obs_cols=()):
    """Load a registered Tabula Sapiens tissue: AnnData with X = ln(1 + CP10k) of the raw counts.

    X is the raw counts (``raw/X``) scaled to 10,000 per cell, where library sizes use all genes, then
    log1p, as float32. ``var_names`` are gene symbols (``raw/var/feature_name``; duplicates keep the
    first), ``obs`` holds ``cell_type``, ``donor_id`` and ``label`` (= cell_type, the groupby), and
    cell types with fewer than ``min_cells`` cells are dropped, as LIANA does.

    ``obs_cols``: further ``obs`` columns to read (e.g. ``broad_cell_class``, ``free_annotation``).

    ``genes``: keep only these genes (e.g. a tool's database). Library sizes still use all genes, so
    the values are the same as without it. The file is read with h5py, because it also holds a
    ``scale_data`` layer too large for 8 GB of memory.
    """
    import anndata as ad
    import h5py
    import numpy as np
    import pandas as pd
    from scipy import sparse

    if not name.startswith("tabula_sapiens_"):
        raise ValueError(f"{name} is not a Tabula Sapiens tissue; registered: "
                         f"{sorted(n for n in DATASETS if n.startswith('tabula_sapiens_'))}")
    with h5py.File(path(name), "r") as f:
        g = f["raw/X"]
        n_obs, n_var = len(g["indptr"]) - 1, len(f["raw/var/_index"])
        counts = sparse.csr_matrix((g["data"][:], g["indices"][:].astype(np.int32), g["indptr"][:]),
                                   shape=(n_obs, n_var))
        symbols = _h5_strings(f, "raw/var/feature_name")
        obs = pd.DataFrame({"cell_type": _h5_strings(f, "obs/cell_type"),
                            "donor_id": _h5_strings(f, "obs/donor_id"),
                            **{c: _h5_strings(f, f"obs/{c}") for c in obs_cols}}, index=_h5_strings(f, "obs/_index"))
    lib = np.asarray(counts.sum(1, dtype=np.float64)).ravel()
    keep = np.arange(n_var) if genes is None else np.where(pd.Index(symbols).isin(set(genes)))[0]
    X = sparse.csr_matrix(counts[:, keep], dtype=np.float64)
    del counts
    X = (sparse.diags(1e4 / np.where(lib > 0, lib, 1)) @ X).tocsr()
    X.data = np.log1p(X.data)
    adata = ad.AnnData(X=X.astype(np.float32), obs=obs)
    adata.var_names = pd.Index(symbols[keep]).astype(str)
    adata = adata[:, ~adata.var_names.duplicated()].copy()
    adata.obs["label"] = adata.obs["cell_type"].astype(str)
    vc = adata.obs.label.value_counts()
    adata = adata[adata.obs.label.isin(vc[vc >= min_cells].index)].copy()
    adata.obs["label"] = adata.obs.label.astype("category")
    return adata


if __name__ == "__main__":
    for n in sys.argv[1:] or ["tnbc_wu2021"]:
        print(n, "->", download(n))
