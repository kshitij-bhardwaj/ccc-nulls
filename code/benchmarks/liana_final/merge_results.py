"""Merge result files produced on several machines (or Kaggle sessions) into this checkout's.

    python -m benchmarks.liana_final.merge_results <dir> --machine "Kaggle CPU, 32 GB"

``<dir>`` mirrors ``benchmarks/results/liana_final``: the plan's files at its top, AMENDMENT_2's in
``amendment2/``, AMENDMENT_3's in ``amendment3/`` and AMENDMENT_4's in ``amendment4/``. Rules:

- every file must carry the same plan and amendment SHA-256 as this checkout's, else this stops;
- AMENDMENT_2, 3 and 4 files hold one record per run (``--only`` shards): runs are merged and tagged with
  ``machine``; a run present on both sides must be identical;
- the plan's ``calibration.json`` holds one entry per tissue (``--dataset``): tissues are merged;
- any other file is copied if absent here; if present, it must say the same apart from ``written``, wall times
  (``seconds*`` keys, anywhere) and SHA-256s one side lacks; the same holds for runs and tissues; when the file here has no ``machine``, the labelled one replaces it;
- ``.csv`` files are copied if absent here, else must be byte-identical.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

from .amendment2 import OUT2
from .pipeline import OUT

SHAS = ("plan_sha256", "amendment_1_sha256", "amendment_2_sha256", "amendment_3_sha256", "amendment_4_sha256",
        "amendment_5_sha256", "amendment_6_sha256", "amendment_7_sha256", "amendment_8_sha256",
        "amendment_9_sha256")


def _strip(d: dict, drop=("written", "machine", "merged_from")) -> dict:
    return {k: v for k, v in d.items() if k not in drop}


def _noclock(x):
    """x without wall times (keys seconds*), which differ between machines."""
    if isinstance(x, dict):
        return {k: _noclock(v) for k, v in x.items() if not k.startswith("seconds")}
    return [_noclock(v) for v in x] if isinstance(x, list) else x


def _content(d: dict) -> dict:
    """What a result says: without labels, wall times and the SHA-256s (checked separately)."""
    return _noclock({k: v for k, v in _strip(d).items() if k not in SHAS})


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
            if _noclock(_strip(old["runs"][tag])) != _noclock(_strip(rec)):
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
        new = {"machine": machine} | new                       # a label from further upstream wins
        if f.name == "calibration.json":
            new["machine_by_tissue"] = {t: new["machine"] for t in new["results"]} | new.get("machine_by_tissue", {})
        dst.write_text(json.dumps(new, indent=1, default=str))
        return 1
    old = json.loads(dst.read_text())
    _check(f.name, old, new, f.parent)
    if f.name == "calibration.json":
        added = [t for t in new["results"] if t not in old["results"]]
        labelled = old.get("machine_by_tissue", {})
        for t in new["results"]:
            if t in old["results"] and _noclock(old["results"][t]) != _noclock(new["results"][t]):
                raise SystemExit(f"calibration.json: {t} differs between the two sources")
            if t in old["results"] and t not in labelled:          # same result; keep the labelled source
                added.append(t)
        old["results"].update({t: new["results"][t] for t in added})
        if added:
            old.setdefault("machine_by_tissue", {}).update({t: machine for t in added})
        dst.write_text(json.dumps(old, indent=1, default=str))
        return len(added)
    if _content(old) != _content(new):
        raise SystemExit(f"{f.name}: differs between the two sources")
    if "machine" not in old:                # same result, but only this source says where it ran: keep it
        dst.write_text(json.dumps({"machine": machine} | new, indent=1, default=str))
        return 1
    return 0


def merge(src: Path, machine: str) -> None:
    for sub, dst_dir in (("", OUT), ("amendment2", OUT2), ("amendment3", OUT / "amendment3"), ("amendment4", OUT / "amendment4"),
                         ("amendment5", OUT / "amendment5"),
                         ("amendment6", OUT / "amendment6"), ("amendment7", OUT / "amendment7"),
                         ("amendment8", OUT / "amendment8"), ("amendment9", OUT / "amendment9")):
        dst_dir.mkdir(parents=True, exist_ok=True)
        for f in sorted((src / sub).glob("*.json")):
            dst = dst_dir / f.name
            n = merge_runs(f, dst, machine) if "runs" in json.loads(f.read_text()) and sub else merge_file(f, dst, machine)
            print(f"{sub + '/' if sub else ''}{f.name}: {n} merged from {machine}")
        for f in sorted((src / sub).glob("*.csv")):              # A2.9's and AMENDMENT_4's case tables
            if not (dst_dir / f.name).exists():
                shutil.copyfile(f, dst_dir / f.name)
            elif (dst_dir / f.name).read_bytes() != f.read_bytes():
                raise SystemExit(f"{sub}/{f.name}: differs between the two sources")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_final.merge_results")
    ap.add_argument("src", type=Path)
    ap.add_argument("--machine", required=True)
    args = ap.parse_args(argv)
    merge(args.src, args.machine)


if __name__ == "__main__":
    main()
