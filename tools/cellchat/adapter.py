"""CellChat adapter: dataset -> CellChat inputs -> CellChat run -> common schema. Python only.

"Native" here is ``cellchat_py``: CellChat 2.1.2 (jinworks/CellChat, release tag v2.1.2) ported to Python
line by line, with CellChatDB read from the package's own data file and R's random-number
generator reproduced, so its output is what CellChat itself reports. No R is needed.
Contract: tools/README.md. Common schema: ccc_nulls.schema.

Null B is not implemented here: ``side_score`` splits CellChat's score into its sender and receiver
sides, and ``null_b`` hands that to the shared engine (``ccc_nulls.composite``), so any Null B method
registered there runs on CellChat's exact score with no change to this folder.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ccc_nulls import composite, schema
from ccc_nulls import liana_reference as _liana
from tools.cellchat import cellchat_py as cc
from tools.cellchat import nulls

TOOL = "cellchat"
LIANA_METHOD = "cellchat"
SCORE_NAME = "lr_probs"
NATIVE_SCORE_NAME = "prob"
HAS_NATIVE_NULL = True
CELLCHAT_REF = f"jinworks/CellChat@{cc.RELEASE} ({cc.VERSION}), ported to Python"
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
               population_size: bool = False, version: str = cc.VERSION, log=None, **kw) -> dict:
    """``computeCommunProb(type = "triMean")`` with CellChat's defaults (nboot = 100,
    seed.use = 1, Kh = 0.5, n = 1, population.size = FALSE): prob, pval, permutations.
    ``version``: "2.1.2" (release, default), "2.2.0.9001" (main) or "1.0.0" (the 2021 paper's
    version), all on the same inputs."""
    res = cc.compute_commun_prob(inputs["X"], inputs["genes"], inputs["codes"], inputs["levels"],
                                 inputs["LRsig"], inputs["db"], nboot=nboot, seed=seed,
                                 population_size=population_size, log=log, version=version)
    return {**res, "db": inputs["db"], "nboot": nboot, "seed": seed, "version": version}


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


# ------------------------------------------------------------------ Null B through the shared engine
def side_score(adata, interactions: pd.DataFrame, groupby: str = "label", *,
               genes_signaling: list[str] | None = None, version: str = cc.VERSION) -> dict:
    """CellChat's score as a ``ccc_nulls.composite.SideScore``: the hook between CellChat and Null B.

    ``summarize(codes, K)``: per-group trimeans of data.signaling / max (``nulls.TriMean``).
    ``score(sender, receiver)``: CellChat's prob, nLR x K x K, with the ligand and the sender's
    agonist / antagonist factor read from ``sender`` and the receptor, co-receptors and the receiver's
    factor from ``receiver`` (``nulls.CellChatScore``). Equal arguments give CellChat's prob exactly.

    Returns dict(side, codes, levels, names, data_max): ``names`` are the L-R pairs along the score's
    first axis. Version 2.x only: 1.0.0 forms complexes per cell, so its sides do not split.
    """
    if version not in ("2.1.2", "2.2.0.9001"):
        raise NotImplementedError("CellChat's score splits into sides for versions 2.x; 1.0.0 forms complexes per cell")
    levels, codes, K, lr, U, gene_pos, data_max = nulls._setup(adata, interactions, groupby, genes_signaling)
    tm = nulls.TriMean(U, quantile="stats" if version == "2.1.2" else "collapse")
    sc = nulls.CellChatScore(lr, gene_pos)
    sizes = np.bincount(codes, minlength=K).astype(float)
    side = composite.SideScore(summarize=tm, score=lambda s, r: sc(s, sizes, avg_receiver=r), name=f"{TOOL} {version}")
    return dict(side=side, codes=codes, levels=levels, names=sc.names, data_max=data_max)


def null_b(adata, interactions: pd.DataFrame, groupby: str = "label", *, method="composite", B: int = 1000,
           seed: int = 0, rng: str = "numpy", permutations: np.ndarray | None = None,
           genes_signaling: list[str] | None = None, delta: float | tuple = 0.0, version: str = cc.VERSION,
           **options) -> pd.DataFrame:
    """Null B on CellChat's score: ``ccc_nulls.composite.get(method)`` applied to ``side_score``.

    ``method``: a Null B registered in ``ccc_nulls.composite`` ("composite", the intersection-union
    test, by default) or a callable with the same signature. ``B, seed, rng, permutations`` as
    ``nulls.null_a`` (equal values give the same permutations). Ties are counted as CellChat counts
    them (prob_b > prob, ``ties="ignore"``); ``options`` go to the method.

    Returns ``interactions`` plus the method's columns, read at each row's (pair, source, target):
    for "composite", ``score_null_b`` (CellChat's prob), ``n_pos_{L,R}_null_b``, ``p_L_null_b``,
    ``p_R_null_b``, ``p_null_b`` (first delta), ``limiting_side_null_b`` and ``p_null_b_delta=<d>``.
    Correct with BH over ``testable(null_a_result, null_b_result)``.
    """
    h = side_score(adata, interactions, groupby, genes_signaling=genes_signaling, version=version)
    if permutations is None and rng == "cellchat":
        permutations = cc.r_permutations(len(h["codes"]), B, seed)
    if permutations is not None:
        permutations = np.asarray(permutations) - 1         # CellChat's permutations are 1-based (R)
    options.setdefault("ties", "ignore")
    res = composite.get(method)(h["side"], h["codes"], len(h["levels"]), B=B, seed=seed,
                                permutations=permutations, delta=delta, **options)
    out = interactions.copy()
    for k, v in composite.to_columns(res, nulls._rows(out, h["names"], h["levels"])).items():
        out[k] = v
    out.attrs["null_b"] = dict(method=method if isinstance(method, str) else getattr(method, "__name__", "callable"),
                               test="intersection-union, p = max(p_L, p_R)" if method == "composite" else None,
                               version=version, B=res.get("B", B), seed=int(seed),
                               rng=rng if permutations is None or rng == "cellchat" else "given",
                               deltas=res.get("deltas"), data_max=h["data_max"], groups=[str(x) for x in h["levels"]])
    return out


def testable(null_a_result: pd.DataFrame, null_b_result: pd.DataFrame) -> np.ndarray:
    """Rows to correct over (``ccc_nulls.composite.testable``): some labelling, observed or permuted
    for Null A or either side for Null B, gives prob > 0. Outside, p = 1 under both nulls."""
    b = null_b_result
    return composite.testable(b.score_null_b, b.n_pos_L_null_b, b.n_pos_R_null_b,
                              null_a_result.perm_npos_null_a.to_numpy() > 0)
