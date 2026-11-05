"""Merge result files produced on several machines (or Kaggle sessions) into this checkout's.

    python -m benchmarks.liana_final.merge_results <dir> --machine "Kaggle CPU, 32 GB"

``<dir>`` mirrors ``benchmarks/results/liana_final``: the plan's files at its top, AMENDMENT_2's in
``amendment2/``. Rules:

- every file must carry the same plan and amendment SHA-256 as this checkout's, else this stops;
- AMENDMENT_2 files hold one record per run (``--only`` shards): runs are merged and tagged with
  ``machine``; a run present on both sides must be identical;
- the plan's ``calibration.json`` holds one entry per tissue (``--dataset``): tissues are merged;
- any other file is copied if absent here; if present, it must be identical apart from ``written``.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from .amendment2 import OUT2
from .pipeline import OUT

SHAS = ("plan_sha256", "amendment_1_sha256", "amendment_2_sha256")


def _strip(d: dict, drop=("written", "machine", "merged_from")) -> dict:
    return {k: v for k, v in d.items() if k not in drop}


def _check(name: str, old: dict, new: dict, src: Path):
    for k in SHAS:
        if k in old and k in new and old[k] != new[k]:
            raise SystemExit(f"{name}: {k} differs ({old[k]} here, {new[k]} in {src})")


def merge_runs(f: Path, dst: Path, machine: str) -> int:
    new = json.loads(f.read_text())
    old = json.loads(dst.read_text()) if dst.exists() else {k: new[k] for k in SHAS if k in new} | {"runs": {}}
    _check(f.name, old, new, f.parent)
    added = 0
    for tag, rec in new["runs"].items():
        rec = rec | {"machine": rec.get("machine", machine)}
        if tag in old["runs"]:
            if _strip(old["runs"][tag]) != _strip(rec):
                raise SystemExit(f"{f.name}: run {tag} differs between the two sources")
            continue
        old["runs"][tag] = rec
        added += 1
    for k, v in new.items():
        if k not in old and k != "runs":
            old[k] = v
    if added:
        old.setdefault("merged_from", []).append({"machine": machine, "runs": added, "written": new.get("written")})
    dst.write_text(json.dumps(old, indent=1, default=str))
    return added


def merge_file(f: Path, dst: Path, machine: str) -> int:
    new = json.loads(f.read_text())
    if not dst.exists():
        dst.write_text(json.dumps(new | {"machine": machine}, indent=1, default=str))
        return 1
    old = json.loads(dst.read_text())
    _check(f.name, old, new, f.parent)
    if f.name == "calibration.json":
        added = [t for t in new["results"] if t not in old["results"]]
        for t in new["results"]:
            if t in old["results"] and old["results"][t] != new["results"][t]:
                raise SystemExit(f"calibration.json: {t} differs between the two sources")
        old["results"].update({t: new["results"][t] for t in added})
        if added:
            old.setdefault("machine_by_tissue", {}).update({t: machine for t in added})
        dst.write_text(json.dumps(old, indent=1, default=str))
        return len(added)
    if _strip(old) != _strip(new):
        raise SystemExit(f"{f.name}: differs between the two sources")
    return 0


def merge(src: Path, machine: str) -> None:
    for sub, dst_dir in (("", OUT), ("amendment2", OUT2)):
        dst_dir.mkdir(parents=True, exist_ok=True)
        for f in sorted((src / sub).glob("*.json")):
            dst = dst_dir / f.name
            n = merge_runs(f, dst, machine) if "runs" in json.loads(f.read_text()) and sub else merge_file(f, dst, machine)
            print(f"{sub + '/' if sub else ''}{f.name}: {n} merged from {machine}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.merge_results")
    ap.add_argument("src", type=Path)
    ap.add_argument("--machine", required=True)
    args = ap.parse_args(argv)
    merge(args.src, args.machine)


if __name__ == "__main__":
    main()
