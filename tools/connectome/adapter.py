"""Connectome adapter: dataset -> native inputs -> native run -> common schema.

Contract: tools/README.md. Common schema: ccc_nulls.schema.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ccc_nulls import liana_reference as _liana
from ccc_nulls import schema

TOOL = "connectome"
LIANA_METHOD = "connectome"
SCORE_NAME = "scaled_weight"
HAS_NATIVE_NULL = False


def liana_reference(adata, groupby: str = "label", **kw) -> pd.DataFrame:
    """LIANA's re-implementation of Connectome, in the common schema (working baseline)."""
    return _liana.run(adata, LIANA_METHOD, groupby, **kw)


def prepare_inputs(adata, groupby: str, out_dir: Path) -> dict:
    """Write the files Connectome expects and return their paths.

    Expected inputs (verify against the installed version):
    1. Expression for a Seurat object: counts or log-normalised matrix + cell identities.
    2. Seurat preprocessing: `NormalizeData`, `ScaleData` (on the LR genes at least), `Idents()` = cell type.
    """
    raise NotImplementedError("connectome: prepare_inputs is not implemented yet (see tools/connectome/README.md)")


def run_native(inputs: dict, out_dir: Path, **kw) -> dict:
    """Run the original Connectome on ``inputs`` and return paths or tables of its output.

    How it is run (verify):
    `CreateConnectome(object, species = "human", p.values = TRUE, calculate.DOR = FALSE)`; write the edge table to CSV.
    """
    raise NotImplementedError("connectome: run_native is not implemented yet (see tools/connectome/README.md)")


def native_to_schema(native: dict) -> pd.DataFrame:
    """Convert Connectome's output into the common schema (``score`` = ``scaled_weight``-equivalent).

    Native output (verify):
    Edge table with `source, target, ligand, receptor, weight_norm, weight_sc, p_val_lig, p_val_rec, ...`.

    Must end with ``return schema.validate(df)``.
    """
    raise NotImplementedError("connectome: native_to_schema is not implemented yet (see tools/connectome/README.md)")
