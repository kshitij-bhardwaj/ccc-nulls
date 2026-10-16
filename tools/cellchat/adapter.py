"""CellChat adapter: dataset -> CellChat inputs -> CellChat run -> common schema. Python only.

"Native" here is ``cellchat_py``: CellChat 2.2.0.9001 (jinworks/CellChat@75253cd) ported to Python
line by line, with CellChatDB read from the package's own data file and R's random-number
generator reproduced, so its output is what CellChat itself reports. No R is needed.
Contract: tools/README.md. Common schema: ccc_nulls.schema.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccc_nulls import liana_reference as _liana
from ccc_nulls import schema
from tools.cellchat import cellchat_py as cc

TOOL = "cellchat"
LIANA_METHOD = "cellchat"
SCORE_NAME = "lr_probs"
NATIVE_SCORE_NAME = "prob"
HAS_NATIVE_NULL = True
CELLCHAT_REF = f"jinworks/CellChat@{cc.COMMIT} ({cc.VERSION}), ported to Python"
SEP = ";"


def liana_reference(adata, groupby: str = "label", **kw) -> pd.DataFrame:
    """LIANA's re-implementation of CellChat, in the common schema (working baseline)."""
    return _liana.run(adata, LIANA_METHOD, groupby, **kw)


def prepare_inputs(adata, groupby: str, out_dir: Path | None = None) -> dict:
    """What ``computeCommunProb`` starts from, built as CellChat's tutorial pipeline does:

    ``subsetDB(CellChatDB.human)`` -> ``subsetData`` (data.signaling = the CellChatDB genes of
    ``adata.X``, taken as the normalised data) -> ``identifyOverExpressedGenes`` ->
    ``identifyOverExpressedInteractions`` (LRsig). Writes LRsig.csv and the gene lists to
    ``out_dir`` if given.
    """
    db = cc.subset_db(cc.load_db("human"))
    X, genes = cc.data_signaling(adata, db)
    codes, levels = cc.labels_and_levels(adata.obs[groupby])
    features = [genes[j] for j in cc.identify_overexpressed_genes(X, codes, len(levels))]
    lr = cc.identify_overexpressed_interactions(db, genes, features)
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        lr.to_csv(out_dir / "LRsig.csv")
        (out_dir / "genes_signaling.txt").write_text("\n".join(genes) + "\n")
        (out_dir / "features_sig.txt").write_text("\n".join(features) + "\n")
    return dict(X=X, genes=genes, codes=codes, levels=levels, db=db, LRsig=lr, features_sig=features)


def run_native(inputs: dict, out_dir: Path | None = None, *, nboot: int = 100, seed: int = 1,
               population_size: bool = False, log=None, **kw) -> dict:
    """``computeCommunProb(type = "triMean")`` with CellChat's defaults (nboot = 100,
    seed.use = 1, Kh = 0.5, n = 1, population.size = FALSE): prob, pval, permutations."""
    res = cc.compute_commun_prob(inputs["X"], inputs["genes"], inputs["codes"], inputs["levels"],
                                 inputs["LRsig"], inputs["db"], nboot=nboot, seed=seed,
                                 population_size=population_size, log=log)
    return {**res, "db": inputs["db"], "nboot": nboot, "seed": seed}


def _gene_lists(lr: pd.DataFrame, db: dict, genes_signaling) -> pd.DataFrame:
    """CellChat's ligand / receptor / cofactor identifiers resolved to gene lists."""
    gene_set = set(genes_signaling)
    out = pd.DataFrame(index=lr.index)
    unit = lambda x: [x] if x in gene_set else cc._subunits(db["complex"], x)
    out["ligand_genes"] = [SEP.join(unit(x)) for x in lr.ligand]
    out["receptor_genes"] = [SEP.join(unit(x)) for x in lr.receptor]
    for col, name in [("agonist", "agonist_genes"), ("antagonist", "antagonist_genes"),
                      ("co_A_receptor", "co_A_genes"), ("co_I_receptor", "co_I_genes")]:
        out[name] = [SEP.join(cc._cofactors(db["cofactor"], x)) if x else "" for x in lr[col]]
    return out


def native_to_schema(native: dict) -> pd.DataFrame:
    """``net$prob`` / ``net$pval`` (unfiltered) as the common schema: one row per
    (L-R pair, source, target), including prob = 0 rows (p = 1).

    ``ligand`` / ``receptor`` are CellChat's identifiers (a complex keeps its CellChat name);
    ``ligand_complex`` / ``receptor_complex`` list its genes joined by "_". Extra columns:
    ``interaction_name, pathway_name, annotation`` and the gene lists ``nulls.null_a`` uses.
    """
    lr, groups, prob, pval = native["LRsig"], native["groups"], native["prob"], native["pval"]
    genes = _gene_lists(lr, native["db"], native["genes_signaling"])
    nLR, K = len(lr), len(groups)
    i, s, t = (a.ravel() for a in np.meshgrid(np.arange(nLR), np.arange(K), np.arange(K), indexing="ij"))
    g = np.asarray(groups, dtype=object)
    df = pd.DataFrame({
        "ligand": lr.ligand.to_numpy()[i], "receptor": lr.receptor.to_numpy()[i],
        "source": g[s], "target": g[t],
        "score": prob[i, s, t], "score_name": NATIVE_SCORE_NAME, "p_native": pval[i, s, t],
        "ligand_complex": genes.ligand_genes.str.replace(SEP, "_").to_numpy()[i],
        "receptor_complex": genes.receptor_genes.str.replace(SEP, "_").to_numpy()[i],
        "tool": TOOL,
        "interaction_name": lr.index.to_numpy()[i],
        "pathway_name": lr.pathway_name.to_numpy()[i], "annotation": lr.annotation.to_numpy()[i],
    })
    for c in genes.columns:
        df[c] = genes[c].to_numpy()[i]
    return schema.validate(df, require_p_native=True)
