"""End-to-end with liana-py on scanpy's bundled PBMC subset (no download).

pip install "ccc-nulls[liana]"
python examples/liana_pbmc.py

Note: pbmc68k_reduced keeps only 765 highly variable genes, so the decoy strata hold a
handful of genes and Null B cannot reach small enough p-values to survive BH. The
example shows the wiring and the resolution diagnostic, not a biological result. Run on
the full transcriptome (e.g. a Tabula Sapiens tissue) for real analyses.
"""
import warnings

import liana as li
import scanpy as sc

import ccc_nulls as cn

warnings.filterwarnings("ignore", category=FutureWarning)

adata = sc.datasets.pbmc68k_reduced()
li.mt.cellphonedb(adata, groupby="bulk_labels", expr_prop=0.1, use_raw=True,
                  n_perms=1000, key_added="cpdb", verbose=False)

cands = cn.adapters.from_liana(adata.uns["cpdb"], method="cellphonedb")
resource = li.resource.select_resource("consensus")
resource_genes = set(resource.ligand.str.split("_").explode()) | set(resource.receptor.str.split("_").explode())

res = cn.significance_layer_anndata(
    adata, "bulk_labels", cands, use_raw=True,          # the matrix liana scored
    score="mean", B=1000, resource_genes=resource_genes,
)

print(f"{len(res):,} interactions ({adata.raw.n_vars} genes, {adata.obs.bulk_labels.nunique()} cell types)")
print(cn.overlap_summary(res).to_string())
print(f"\nNull B resolution: median smallest attainable p = {res.p_strat_min.median():.1e}; "
      f"BH needs <= {0.05 / len(res):.1e} for the top hit")
cols = ["ligand_complex", "receptor_complex", "source", "target", "q_rand", "q_strat"]
print("\nCalled under both nulls (top 10 by q_strat):")
print(res[res.flag == "both"].sort_values("q_strat")[cols].head(10).to_string(index=False))
print("\nA-only (abundance-driven), top 10 by q_rand:")
print(res[res.flag == "A_only"].sort_values("q_rand")[cols].head(10).to_string(index=False))
