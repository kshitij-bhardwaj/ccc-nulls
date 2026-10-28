"""
Connectome Adapter for ccc-nulls repository (tools/connectome/adapter.py)
Follows the exact interface of tools/cellchat/adapter.py for canonical schema integration.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad
import liana as li

SEP = ";"

def prepare_inputs(adata: ad.AnnData, groupby: str = "cell_type") -> dict:
    """Prepares AnnData inputs for Connectome."""
    if "normalized" in adata.layers:
        X = adata.layers["normalized"]
    elif sparse.issparse(adata.X):
        X = adata.X
    else:
        X = sparse.csr_matrix(adata.X)
    return {"adata": adata, "groupby": groupby, "X": X}

def run_native(
    inputs: dict,
    resource_name: str = "consensus",
    expr_prop: float = 0.1,
    min_cells: int = 5,
    **kwargs
) -> pd.DataFrame:
    """
    Executes Connectome via LIANA.
    Returns native LIANA Connectome DataFrame.
    """
    adata = inputs["adata"]
    groupby = inputs["groupby"]
    
    # Run LIANA Connectome
    df = li.method.connectome(
        adata=adata,
        groupby=groupby,
        resource_name=resource_name,
        expr_prop=expr_prop,
        min_cells=min_cells,
        n_perms=1,         # Permutations handled externally by Null A & Null B engines
        use_raw=False,
        verbose=False
    )
    return df

def native_to_schema(native_df: pd.DataFrame) -> pd.DataFrame:
    """
    Maps native LIANA Connectome DataFrame into the shared common-schema:
    Columns:
      interaction_name, source, target, ligand, receptor,
      ligand_complex, receptor_complex, ligand_genes, receptor_genes,
      score (expr_prod)
    """
    if native_df is None or len(native_df) == 0:
        return pd.DataFrame()
        
    df = native_df.copy()
    
    # Ensure source and target
    src_col = "source" if "source" in df.columns else df.columns[0]
    tgt_col = "target" if "target" in df.columns else df.columns[1]
    
    # Ligand / Receptor
    lig_col = "ligand_complex" if "ligand_complex" in df.columns else "ligand"
    rec_col = "receptor_complex" if "receptor_complex" in df.columns else "receptor"
    
    # Score
    score_col = "expr_prod" if "expr_prod" in df.columns else "scaled_weight"
    
    df["interaction_name"] = df[lig_col].astype(str) + "|" + df[rec_col].astype(str)
    df["source"] = df[src_col].astype(str)
    df["target"] = df[tgt_col].astype(str)
    df["ligand"] = df[lig_col].astype(str)
    df["receptor"] = df[rec_col].astype(str)
    df["ligand_complex"] = df[lig_col].astype(str)
    df["receptor_complex"] = df[rec_col].astype(str)
    df["ligand_genes"] = df[lig_col].astype(str).str.replace("_", SEP)
    df["receptor_genes"] = df[rec_col].astype(str).str.replace("_", SEP)
    df["score"] = df[score_col].astype(float)
    
    cols = [
        "interaction_name", "source", "target", "ligand", "receptor",
        "ligand_complex", "receptor_complex", "ligand_genes", "receptor_genes", "score"
    ]
    return df[[c for c in cols if c in df.columns]]
