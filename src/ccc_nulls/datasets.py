"""Dataset registry, downloads and loaders.

Files live in ``<repo>/data/`` by default; set ``CCC_NULLS_DATA`` to use another folder.

    python -m ccc_nulls.datasets tnbc_wu2021          # download + md5 check
    >>> from ccc_nulls.datasets import load_openproblems
    >>> adata, truth = load_openproblems("tnbc_wu2021")
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
        url="https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5",
        filename=None,
        md5=None,
        size=None,
        organism="human",
        groupby="cell_type",
        sample_key="donor_id",
        truth_keys=None,
        summary="Tabula Sapiens (CELLxGENE). Download one tissue h5ad from the collection page, save it "
                "as data/tabula_sapiens_<tissue>.h5ad, and pass that path to load_h5ad().",
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


if __name__ == "__main__":
    for n in sys.argv[1:] or ["tnbc_wu2021"]:
        print(n, "->", download(n))
