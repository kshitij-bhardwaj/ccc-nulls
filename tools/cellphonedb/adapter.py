"""CellPhoneDB (v5) adapter: dataset -> native inputs -> native run -> common schema.

Contract: tools/README.md. Common schema: ccc_nulls.schema.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ccc_nulls import liana_reference as _liana
from ccc_nulls import schema

TOOL = "cellphonedb"
LIANA_METHOD = "cellphonedb"
SCORE_NAME = "lr_means"
HAS_NATIVE_NULL = True


def liana_reference(adata, groupby: str = "label", **kw) -> pd.DataFrame:
    """LIANA's re-implementation of CellPhoneDB (v5), in the common schema (working baseline)."""
    return _liana.run(adata, LIANA_METHOD, groupby, **kw)


def prepare_inputs(adata, groupby: str, out_dir: Path) -> dict:
    """Write the files CellPhoneDB (v5) expects and return their paths.

    Expected inputs (verify against the installed version):
    1. `counts.h5ad`: log-normalised expression, cells x genes, gene symbols as var_names.
    2. `meta.tsv`: two columns `barcode_sample`, `cell_type` (cell barcode -> label).
    3. The CellPhoneDB database zip (`cellphonedb.zip`), downloaded with
       `cellphonedb.utils.db_utils.download_database(out_dir, "v5.0.0")`.
    """
    raise NotImplementedError("cellphonedb: prepare_inputs is not implemented yet (see tools/cellphonedb/README.md)")


def run_native(inputs: dict, out_dir: Path, **kw) -> dict:
    """Run the original CellPhoneDB (v5) on ``inputs`` and return paths or tables of its output.

    How it is run (verify):
    `cpdb_statistical_analysis_method.call(cpdb_file_path=..., meta_file_path=..., counts_file_path=...,
    counts_data="hgnc_symbol", output_path=..., iterations=1000, threshold=0.1)`
    """
    raise NotImplementedError("cellphonedb: run_native is not implemented yet (see tools/cellphonedb/README.md)")


def native_to_schema(native: dict) -> pd.DataFrame:
    """Convert CellPhoneDB (v5)'s output into the common schema (``score`` = ``lr_means``-equivalent).

    Native output (verify):
    Wide tables `means`, `pvalues`, `significant_means`: one row per interacting pair
    (`interacting_pair`, `partner_a`, `partner_b`, `gene_a`, `gene_b`, `receptor_a`, `receptor_b`, ...)
    and one column per `"<source>|<target>"`. Partner order is **not** ligand-first: use the
    `receptor_a` / `receptor_b` flags to decide which partner is the ligand.

    Must end with ``return schema.validate(df)``.
    """
    raise NotImplementedError("cellphonedb: native_to_schema is not implemented yet (see tools/cellphonedb/README.md)")
