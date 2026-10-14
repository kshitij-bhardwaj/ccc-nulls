"""CellChat adapter: dataset -> native inputs -> native run -> common schema.

Contract: tools/README.md. Common schema: ccc_nulls.schema.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ccc_nulls import liana_reference as _liana
from ccc_nulls import schema

TOOL = "cellchat"
LIANA_METHOD = "cellchat"
SCORE_NAME = "lr_probs"
HAS_NATIVE_NULL = True


def liana_reference(adata, groupby: str = "label", **kw) -> pd.DataFrame:
    """LIANA's re-implementation of CellChat, in the common schema (working baseline)."""
    return _liana.run(adata, LIANA_METHOD, groupby, **kw)


def prepare_inputs(adata, groupby: str, out_dir: Path) -> dict:
    """Write the files CellChat expects and return their paths.

    Expected inputs (verify against the installed version):
    1. Expression: log-normalised genes x cells, written as `matrix.mtx` + `genes.tsv` + `barcodes.tsv`.
    2. `meta.csv`: cell barcode -> `labels`.
    3. Database: `CellChatDB.human` (ships with the R package).
    """
    raise NotImplementedError("cellchat: prepare_inputs is not implemented yet (see tools/cellchat/README.md)")


def run_native(inputs: dict, out_dir: Path, **kw) -> dict:
    """Run the original CellChat on ``inputs`` and return paths or tables of its output.

    How it is run (verify):
    An R script (e.g. `run_cellchat.R`): `createCellChat(object = data, meta = meta, group.by = "labels")`,
    `subsetData`, `identifyOverExpressedGenes`, `identifyOverExpressedInteractions`,
    `computeCommunProb(type = "triMean")`, `filterCommunication(min.cells = 10)`, then
    `subsetCommunication(thresh = 1)` to keep non-significant rows too; write it to CSV.
    """
    raise NotImplementedError("cellchat: run_native is not implemented yet (see tools/cellchat/README.md)")


def native_to_schema(native: dict) -> pd.DataFrame:
    """Convert CellChat's output into the common schema (``score`` = ``lr_probs``-equivalent).

    Native output (verify):
    Long table: `source, target, ligand, receptor, prob, pval, interaction_name, ...`.

    Must end with ``return schema.validate(df)``.
    """
    raise NotImplementedError("cellchat: native_to_schema is not implemented yet (see tools/cellchat/README.md)")
