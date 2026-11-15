"""Paths, logging and result files shared by the parts."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

from ccc_nulls import datasets

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "results" / "liana_nulls"
CORE_DATASETS = ["tnbc_wu2021", "allen_brain_atlas", "tabula_sapiens_mammary", "tabula_sapiens_liver"]
AMENDMENT_DATASETS = ["moffitt2018_hypothalamus"]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def rows_dir() -> Path:
    d = Path(datasets.data_dir()) / "liana_nulls"
    d.mkdir(parents=True, exist_ok=True)
    return d


def cpdb_database() -> Path:
    return Path(datasets.data_dir()) / "cellphonedb_v5.0.0" / "cellphonedb.zip"


def plan_hash() -> str:
    return hashlib.sha256((HERE / "PLAN.md").read_bytes()).hexdigest()


def amendment_hash(n: int = 1) -> str:
    return hashlib.sha256((HERE / f"AMENDMENT_{n}.md").read_bytes()).hexdigest()


def save(name: str, obj: dict, amendment: bool = False):
    """``amendment``: the result belongs to Amendment 1 and records its hash too."""
    OUT.mkdir(parents=True, exist_ok=True)
    head = {"plan_sha256": plan_hash()} | ({"amendment_1_sha256": amendment_hash()} if amendment else {})
    obj = {**head, "written": time.strftime("%Y-%m-%d %H:%M"), **obj}
    (OUT / f"{name}.json").write_text(json.dumps(obj, indent=1, default=str))
    log("wrote", OUT / f"{name}.json")
