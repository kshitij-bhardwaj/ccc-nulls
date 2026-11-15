"""CellChat 2.1.2 run in R (``run_cellchat.R``), and how it compares with the Python port.

The port (``cellchat_py``) is what the nulls use; this module checks it against CellChat itself.
``run_r`` writes ``data.signaling`` as raw float64 (no text rounding), runs the tutorial pipeline in
R (createCellChat -> subsetDB -> subsetData -> identifyOverExpressedGenes ->
identifyOverExpressedInteractions -> computeCommunProb) and returns the same dict as
``adapter.run_native``, so ``adapter.native_to_schema`` reads both. ``fidelity`` compares the two.
Needs R with CellChat 2.1.2 (``remotes::install_github("jinworks/CellChat@v2.1.2")``).
"""
from __future__ import annotations

import hashlib
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from scipy import sparse

from tools.cellchat import cellchat_py as cc

SCRIPT = Path(__file__).with_name("run_cellchat.R")


def write_inputs(adata, groupby: str, out_dir: Path, species: str = "human") -> dict:
    db = cc.subset_db(cc.load_db(species))
    X, genes = cc.data_signaling(adata, db)
    codes, levels = cc.labels_and_levels(adata.obs[groupby])
    G = sparse.csc_matrix(X.T)                      # genes x cells, a column per cell
    G.sort_indices()
    out_dir.mkdir(parents=True, exist_ok=True)
    G.indices.astype("<i4").tofile(out_dir / "i.bin")
    G.indptr.astype("<i4").tofile(out_dir / "p.bin")
    G.data.astype("<f8").tofile(out_dir / "x.bin")
    (out_dir / "genes.txt").write_text("\n".join(genes) + "\n")
    (out_dir / "cells.txt").write_text("\n".join(f"c{i}" for i in range(X.shape[0])) + "\n")
    (out_dir / "labels.txt").write_text("\n".join(levels[codes]) + "\n")
    (out_dir / "levels.txt").write_text("\n".join(levels) + "\n")
    (out_dir / "species.txt").write_text(species + "\n")
    h = hashlib.sha256()
    for f in ("i.bin", "p.bin", "x.bin", "genes.txt", "labels.txt", "levels.txt", "species.txt"):
        h.update((out_dir / f).read_bytes())
    return dict(db=db, genes=genes, levels=list(levels), inputs_sha256=h.hexdigest())


def run_r(adata, groupby: str, out_dir: Path | None = None, *, nboot: int = 100, seed: int = 1,
          species: str = "human", rscript: str = "Rscript", log=None) -> dict:
    """CellChat's own ``net$prob`` / ``net$pval`` on ``adata`` (normalised X), as ``adapter.run_native``.
    With ``out_dir``, an earlier R run there on byte-identical inputs and the same nboot and seed is
    reused (``done.txt`` holds its input hash)."""
    tmp = None
    if out_dir is None:
        tmp = tempfile.TemporaryDirectory()
        out_dir = Path(tmp.name)
    out_dir = Path(out_dir)
    meta = write_inputs(adata, groupby, out_dir, species)
    key = f"{meta['inputs_sha256']} nboot={nboot} seed={seed}"
    done = out_dir / "done.txt"
    if not (done.exists() and done.read_text() == key):
        if log:
            log(f"  R: CellChat on {adata.n_obs:,} cells, nboot = {nboot}, seed = {seed}")
        done.unlink(missing_ok=True)
        proc = subprocess.run([rscript, str(SCRIPT), str(out_dir), str(nboot), str(seed)], capture_output=True, text=True)
        if proc.returncode:
            raise RuntimeError(f"run_cellchat.R failed:\n{proc.stderr[-3000:]}")
        done.write_text(key)
    elif log:
        log(f"  R: reusing the CellChat run in {out_dir}")
    names = (out_dir / "LRsig.txt").read_text().split("\n")[:-1]
    K = len(meta["levels"])
    # R's net$prob is sender x receiver x L-R pair; the port's is L-R pair x sender x receiver
    read = lambda f: np.fromfile(out_dir / f, dtype="<f8").reshape((K, K, len(names)), order="F").transpose(2, 0, 1)
    LRsig = meta["db"]["interaction"].loc[names]
    out = dict(prob=read("prob.bin"), pval=read("pval.bin"), groups=meta["levels"], LRsig=LRsig,
               genes_signaling=meta["genes"], db=meta["db"], nboot=nboot, seed=seed, version="R " + cc.VERSION,
               features_sig=(out_dir / "features_sig.txt").read_text().split("\n")[:-1],
               session=(out_dir / "session.txt").read_text().split("\n")[:-1], inputs_sha256=meta["inputs_sha256"])
    if tmp is not None:
        tmp.cleanup()
    return out


def fidelity(r: dict, py: dict, alpha: float = 0.05) -> dict:
    """R CellChat vs the Python port on the same inputs: candidate genes and L-R pairs, then prob,
    p and CellChat's calls (p < alpha, prob > 0) for every (L-R, sender, receiver)."""
    out = dict(features_sig_identical=set(r["features_sig"]) == set(py.get("features_sig", r["features_sig"])),
               n_lr_r=len(r["LRsig"]), n_lr_python=len(py["LRsig"]),
               lr_identical=list(r["LRsig"].index) == list(py["LRsig"].index),
               groups_identical=list(r["groups"]) == list(py["groups"]))
    if not (out["lr_identical"] and out["groups_identical"]):
        common = r["LRsig"].index.intersection(py["LRsig"].index)
        out.update(n_lr_common=len(common), lr_only_r=sorted(set(r["LRsig"].index) - set(common))[:20],
                   lr_only_python=sorted(set(py["LRsig"].index) - set(common))[:20])
        ir, ip = r["LRsig"].index.get_indexer(common), py["LRsig"].index.get_indexer(common)
        pr, vr, pp, vp = r["prob"][ir], r["pval"][ir], py["prob"][ip], py["pval"][ip]
    else:
        pr, vr, pp, vp = r["prob"], r["pval"], py["prob"], py["pval"]
    call_r, call_p = (vr < alpha) & (pr > 0), (vp < alpha) & (pp > 0)
    out.update(n_rows=int(pr.size), prob_identical=int((pr == pp).sum()),
               prob_max_abs_diff=float(np.abs(pr - pp).max()), pval_identical=int((vr == vp).sum()),
               pval_max_abs_diff=float(np.abs(vr - vp).max()),
               calls_r=int(call_r.sum()), calls_python=int(call_p.sum()), calls_identical=bool((call_r == call_p).all()),
               calls_only_r=int((call_r & ~call_p).sum()), calls_only_python=int((call_p & ~call_r).sum()))
    return out
