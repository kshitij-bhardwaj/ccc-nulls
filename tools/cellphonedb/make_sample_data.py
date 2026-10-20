"""Create a tiny deterministic human-expression input for CellPhoneDB development.

Usage
-----
python tools/cellphonedb/make_sample_data.py

The script writes CellPhoneDB's text counts format and matching metadata to
``data/cellphonedb_sample/``.  The files are intentionally git-ignored because
they are generated artifacts; regenerate them whenever needed.
"""
from __future__ import annotations

import csv
from pathlib import Path


GENES = [
    "CXCL12", "CXCR4", "IL6", "IL6R", "TNF", "TNFRSF1A", "TGFB1", "TGFBR1",
    "CD3D", "MS4A1", "LYZ", "EPCAM",
]
CELL_TYPES = ("Fibroblast", "Myeloid", "T_cell", "B_cell")
CELLS = [f"cell{kind}{i}" for kind in ("F", "M", "T", "B") for i in range(1, 7)]


def expression(gene: str, cell: str) -> float:
    """Small normalized-expression matrix with several detectable LR signals."""
    kind = cell[4]
    base = 0.05
    high = {
        ("CXCL12", "F"): 4.0, ("CXCR4", "T"): 3.5,
        ("IL6", "M"): 3.8, ("IL6R", "T"): 2.8,
        ("TNF", "M"): 3.6, ("TNFRSF1A", "F"): 2.7,
        ("TGFB1", "F"): 3.2, ("TGFBR1", "B"): 2.6,
        ("CD3D", "T"): 3.0, ("MS4A1", "B"): 3.0,
        ("LYZ", "M"): 3.0, ("EPCAM", "F"): 2.0,
    }
    # A small deterministic cell-level change prevents identical columns.
    return high.get((gene, kind), base) + 0.01 * int(cell[-1])


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    out = root / "data" / "cellphonedb_sample"
    out.mkdir(parents=True, exist_ok=True)

    with (out / "counts.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene", *CELLS])
        for gene in GENES:
            writer.writerow([gene, *(f"{expression(gene, cell):.2f}" for cell in CELLS)])

    with (out / "meta.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Cell", "cell_type"])
        for cell in CELLS:
            writer.writerow([cell, CELL_TYPES["FMTB".index(cell[4])]])

    print(f"Wrote {len(GENES)} genes and {len(CELLS)} cells to {out}")


if __name__ == "__main__":
    main()
