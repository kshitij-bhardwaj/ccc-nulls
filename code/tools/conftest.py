"""Shared fixtures for the tool workspaces (small data that ships with scanpy, no download)."""
import pytest


@pytest.fixture(scope="session")
def small_adata():
    """700 PBMCs, 10 cell types, log-normalised expression in X (scanpy's pbmc68k_reduced .raw)."""
    sc = pytest.importorskip("scanpy")
    import anndata as ad

    pb = sc.datasets.pbmc68k_reduced()
    adata = ad.AnnData(X=pb.raw.X.copy(), obs=pb.obs[["bulk_labels"]].copy(), var=pb.raw.var[[]].copy())
    adata.obs["label"] = adata.obs["bulk_labels"].astype(str).astype("category")
    return adata
