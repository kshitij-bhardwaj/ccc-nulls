"""NATMI adapter: dataset -> native inputs -> native run -> common schema.

Contract: tools/README.md. Common schema: ccc_nulls.schema.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ccc_nulls import liana_reference as _liana
from ccc_nulls import schema

TOOL = "natmi"
LIANA_METHOD = "natmi"
SCORE_NAME = "spec_weight"
HAS_NATIVE_NULL = False


def liana_reference(adata, groupby: str = "label", **kw) -> pd.DataFrame:
    """LIANA's re-implementation of NATMI, in the common schema (working baseline)."""
    return _liana.run(adata, LIANA_METHOD, groupby, **kw)


def prepare_inputs(adata, groupby: str, out_dir: Path) -> dict:
    """Write the files NATMI expects and return their paths.

    Expected inputs (verify against the installed version):
    1. `em.csv`: log-normalised expression, genes x cells, first column = gene symbol.
    2. `ann.csv`: two columns, cell barcode and cell-type annotation.
    3. Interaction database: `lrc2p` (connectomeDB2020; ships with NATMI).
    """
    raise NotImplementedError("natmi: prepare_inputs is not implemented yet (see tools/natmi/README.md)")


def run_native(inputs: dict, out_dir: Path, **kw) -> dict:
    """Run the original NATMI on ``inputs`` and return paths or tables of its output.

    How it is run (verify):
    `python ExtractEdges.py --species human --emFile em.csv --annFile ann.csv --interDB lrc2p --out <dir>`
    """
    raise NotImplementedError("natmi: run_native is not implemented yet (see tools/natmi/README.md)")


def native_to_schema(native: dict) -> pd.DataFrame:
    """Convert NATMI's output into the common schema (``score`` = ``spec_weight``-equivalent).

    Native output (verify):
    `Edges_lrc2p.csv`: `Sending cluster, Ligand symbol, Receptor symbol, Target cluster`,
    detection rates, `Edge average expression weight`, `Edge average expression derived specificity`, ...

    Must end with ``return schema.validate(df)``.
    """
    raise NotImplementedError("natmi: native_to_schema is not implemented yet (see tools/natmi/README.md)")
