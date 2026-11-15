"""NATMI adapter: dataset -> NATMI inputs -> NATMI run -> common schema; NATMI's score for the nulls.

"Native" is NATMI itself: ``ExtractEdges.py`` of asrhou/NATMI at ``NATMI_COMMIT``, run in its own
environment (``requirements-native.txt``; the clone is ``$NATMI_DIR``, its interpreter ``$NATMI_PYTHON``).
Contract: tools/README.md. Common schema: ccc_nulls.schema. Plan: benchmarks/liana_nulls/AMENDMENT_5.md.

NATMI's score (``ExtractEdges.py``: ``SplitIntoSinalProteins`` and ``GenSingleCell2CellEdge``), for ligand L
in sender s and receptor R in receiver t, with the means over every cell of a type (zeros included)::

    spec_weight = (mean_L[s] / sum_k mean_L[k]) * (mean_R[t] / sum_k mean_R[k])    "Edge average expression
                                                                                    derived specificity"
    expr_prod   = mean_L[s] * mean_R[t]                                             "Edge average expression weight"

A ratio with a zero denominator is 0 (NATMI's ``fillna(0.0)``). The ligand factor reads only the ligand's
group means and the receptor factor only the receptor's, so the score splits into a sender and a receiver
side exactly (:func:`side_score`), and Null B runs on it through ``ccc_nulls.composite``.

NATMI has no p-value. Null A here is the label permutation the other tools use, applied to NATMI's score,
on the same ``numpy.random.default_rng(seed)`` sequence as Null B. NATMI's database (lrc2p, connectomeDB2020)
has no heteromeric complexes: every ligand and receptor is one gene.
"""
from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

from ccc_nulls import composite, schema
from ccc_nulls import liana_reference as _liana
from ccc_nulls._matrix import as_matrix, cluster_means, columns

TOOL = "natmi"
LIANA_METHOD = "natmi"
LIANA_RESOURCE = "connectomedb2020"          # NATMI's own database, as LIANA+ ships it
SCORE_NAME = "spec_weight"                   # LIANA+'s name for the statistic tested
NATIVE_SCORE_NAME = "Edge average expression derived specificity"
NATIVE_WEIGHT_NAME = "Edge average expression weight"
HAS_NATIVE_NULL = False
NATMI_COMMIT = "f35f677cb2e4c8a176a6c5501a9483c9e28661a1"   # asrhou/NATMI master, 2021-04-13
INTER_DB = "lrc2p"
SCORES = ("spec_weight", "expr_prod")
FIDELITY_RTOL = 1e-9                         # NATMI round-trips its edges through .xlsx files


class NativeNotInstalled(FileNotFoundError):
    """The NATMI clone or its environment is missing (see requirements-native.txt)."""


def natmi_dir() -> Path:
    """The NATMI clone: ``$NATMI_DIR``, else ``<data dir>/natmi/NATMI``."""
    d = os.environ.get("NATMI_DIR")
    if d:
        return Path(d)
    from ccc_nulls import datasets
    return Path(datasets.data_dir()) / "natmi" / "NATMI"


def natmi_python() -> str:
    """The interpreter of NATMI's own environment: ``$NATMI_PYTHON``, else ``python``."""
    return os.environ.get("NATMI_PYTHON", "python")


def load_pairs(path: str | Path | None = None) -> pd.DataFrame:
    """NATMI's lrc2p pairs (``lrdbs/lrc2p.csv`` of the clone) as ``ligand, receptor``, sorted, unique."""
    path = Path(path) if path is not None else natmi_dir() / "lrdbs" / f"{INTER_DB}.csv"
    if not path.exists():
        raise NativeNotInstalled(f"{path} not found: clone NATMI and set NATMI_DIR (tools/natmi/requirements-native.txt)")
    d = pd.read_csv(path, dtype=str)
    d = d.rename(columns={"Ligand gene symbol": "ligand", "Receptor gene symbol": "receptor"})[["ligand", "receptor"]]
    d = d.dropna().apply(lambda c: c.str.strip())
    return d.drop_duplicates().sort_values(["ligand", "receptor"]).reset_index(drop=True)


def database_genes(pairs: pd.DataFrame | None = None) -> set:
    """Every gene of NATMI's database."""
    pairs = load_pairs() if pairs is None else pairs
    return set(pairs.ligand) | set(pairs.receptor)


def side_genes(pairs: pd.DataFrame | None = None) -> tuple[set, set]:
    """Genes that are only ever ligands, and only ever receptors, in NATMI's database (calibration)."""
    pairs = load_pairs() if pairs is None else pairs
    lig, rec = set(pairs.ligand), set(pairs.receptor)
    return lig - rec, rec - lig


def liana_reference(adata, groupby: str = "label", **kw) -> pd.DataFrame:
    """LIANA's re-implementation of NATMI, in the common schema (``spec_weight``), on NATMI's own database."""
    kw.setdefault("resource_name", LIANA_RESOURCE)
    return _liana.run(adata, LIANA_METHOD, groupby, **kw)


# ------------------------------------------------------------------ NATMI's score, split into sides
def labels_and_levels(labels) -> tuple[np.ndarray, list]:
    """Group code per cell and the groups, sorted as NATMI sorts its clusters."""
    lab = np.asarray(labels).astype(str)
    levels = np.unique(lab)
    return np.searchsorted(levels, lab).astype(np.int64), [str(x) for x in levels]


class NatmiScore:
    """``summarize(codes, K)``: types x genes means of the pairs' genes. ``score(sender, receiver)``: the
    score of every (pair, source, target), shape (pairs, K, K), the ligand side read from ``sender`` only
    and the receptor side from ``receiver`` only."""

    def __init__(self, X, genes, pairs: pd.DataFrame, score: str = SCORE_NAME):
        if score not in SCORES:
            raise ValueError(f"score must be one of {SCORES}")
        genes = pd.Index(np.asarray(genes).astype(str))
        self.used = pd.Index(sorted(set(pairs.ligand) | set(pairs.receptor)))
        idx = genes.get_indexer(self.used)
        if (idx < 0).any():
            raise ValueError(f"{int((idx < 0).sum())} genes of the pairs are not in the data")
        self.Xu = columns(as_matrix(X), idx)
        self.li = self.used.get_indexer(pairs.ligand)
        self.ri = self.used.get_indexer(pairs.receptor)
        self.score_name = score

    def summarize(self, codes, n_types: int) -> np.ndarray:
        return cluster_means(self.Xu, np.asarray(codes), n_types)

    def side(self, M: np.ndarray, idx: np.ndarray) -> np.ndarray:
        """types x pairs: the ligand (or receptor) factor of the score."""
        m = M[:, idx]
        if self.score_name == "expr_prod":
            return m
        tot = m.sum(0, keepdims=True)
        return np.divide(m, tot, out=np.zeros_like(m), where=tot > 0)

    def __call__(self, sender: np.ndarray, receiver: np.ndarray) -> np.ndarray:
        a, b = self.side(sender, self.li), self.side(receiver, self.ri)
        return a.T[:, :, None] * b.T[:, None, :]


def side_score(adata, groupby: str = "label", *, pairs: pd.DataFrame | None = None,
               score: str = SCORE_NAME) -> dict:
    """NATMI's score as a ``ccc_nulls.composite.SideScore``: the hook between NATMI and Null B.

    ``pairs``: ``ligand, receptor`` (default: lrc2p), kept where both genes are in the data. Returns
    dict(side, codes, levels, pairs); the score's axes are (pairs, source, target).
    """
    pairs = load_pairs() if pairs is None else pairs
    present = set(adata.var_names.astype(str))
    pairs = pairs[pairs.ligand.isin(present) & pairs.receptor.isin(present)].reset_index(drop=True)
    if len(pairs) == 0:
        raise ValueError("no NATMI pair has both genes in the data")
    codes, levels = labels_and_levels(adata.obs[groupby])
    sc = NatmiScore(adata.X, adata.var_names, pairs, score)
    side = composite.SideScore(summarize=sc.summarize, score=sc, name=f"{TOOL} {score}")
    return dict(side=side, codes=codes, levels=levels, pairs=pairs, score=score)


def grid(h: dict) -> pd.DataFrame:
    """One common-schema row per (pair, source, target), in the order of the score array's ``ravel``."""
    pairs, levels = h["pairs"], np.asarray(h["levels"], dtype=object)
    K = len(levels)
    i, s, t = (a.ravel() for a in np.meshgrid(np.arange(len(pairs)), np.arange(K), np.arange(K), indexing="ij"))
    lig, rec = pairs.ligand.to_numpy(dtype=object)[i], pairs.receptor.to_numpy(dtype=object)[i]
    return pd.DataFrame({"ligand": lig, "receptor": rec, "source": levels[s], "target": levels[t],
                         "ligand_complex": lig, "receptor_complex": rec, "tool": TOOL})


# ------------------------------------------------------------------ the nulls
def null_a(h: dict, *, B: int = 1000, seed: int = 0) -> dict:
    """Null A on NATMI's score: every label permuted, p = (1 + #{s_b >= s}) / (B + 1), 1 where s = 0.
    The permutations are ``default_rng(seed).permutation(codes)``, one per iteration: Null B's sequence
    for the same ``B`` and ``seed`` (``composite.composite``)."""
    side, codes, K = h["side"], np.asarray(h["codes"]), len(h["levels"])
    obs = side.summarize(codes, K)
    s_obs = np.asarray(side.score(obs, obs))
    k = np.zeros(s_obs.shape, np.int64)
    npos = np.zeros(s_obs.shape, np.int64)
    gen = np.random.default_rng(seed)
    for _ in range(B):
        m = side.summarize(gen.permutation(codes), K)
        s = np.asarray(side.score(m, m))
        k += s >= s_obs
        npos += s > 0
    p = np.where(s_obs == 0, 1.0, (1 + k) / (B + 1))
    return dict(score=s_obs, p=p, n_exceed=k, n_pos=npos, B=int(B), seed=int(seed))


def null_b(h: dict, *, method: str = "composite", B: int = 1000, seed: int = 0, delta=0.0, **options) -> dict:
    """Null B on NATMI's score: ``ccc_nulls.composite.get(method)`` on :func:`side_score`. Ties are counted
    (the repository's convention; NATMI has none of its own)."""
    options.setdefault("ties", "count")
    return composite.get(method)(h["side"], h["codes"], len(h["levels"]), B=B, seed=seed, delta=delta, **options)


def testable(a: dict, b: dict) -> np.ndarray:
    """Rows to correct over: some labelling (observed, Null A's, either side of Null B's) gives score > 0."""
    return composite.testable(b["score"], b["n_pos_L"], b["n_pos_R"], a["n_pos"] > 0)


def run(adata, groupby: str = "label", *, B: int = 1000, seed: int = 0, pairs: pd.DataFrame | None = None,
        score: str = SCORE_NAME) -> tuple[pd.DataFrame, dict]:
    """NATMI's score recomputed in Python on every (pair, source, target), with Null A and Null B on the
    same permutations. Columns: the schema's, ``expr_prod``, ``p_null_a``, ``n_exceed_null_a``,
    ``n_pos_null_a``, Null B's (``composite.to_columns``) and ``family``."""
    t0 = time.time()
    h = side_score(adata, groupby, pairs=pairs, score=score)
    df = grid(h)
    a = null_a(h, B=B, seed=seed)
    t1 = time.time()
    b = null_b(h, B=B, seed=seed)
    t2 = time.time()
    if not np.array_equal(a["score"], b["score"]):
        raise RuntimeError("Null B scored something other than Null A")
    df["score"], df["score_name"], df["p_native"] = a["score"].ravel(), score, np.nan
    other = "expr_prod" if score == SCORE_NAME else SCORE_NAME
    sc2 = NatmiScore(adata.X, adata.var_names, h["pairs"], other)
    m = sc2.summarize(h["codes"], len(h["levels"]))
    df[other] = sc2(m, m).ravel()
    df["p_null_a"], df["n_exceed_null_a"], df["n_pos_null_a"] = a["p"].ravel(), a["n_exceed"].ravel(), a["n_pos"].ravel()
    for k, v in composite.to_columns(b).items():
        df[k] = np.asarray(v).ravel()
    df["family"] = testable(a, b).ravel()
    meta = dict(score=score, B=int(B), seed=int(seed), n_pairs=int(len(h["pairs"])), n_types=len(h["levels"]),
                n_rows=int(len(df)), n_family=int(df.family.sum()),
                seconds_null_a=round(t1 - t0, 1), seconds_null_b=round(t2 - t1, 1))
    return df, meta


# ------------------------------------------------------------------ native NATMI
def prepare_inputs(adata, groupby: str, out_dir: Path, *, pairs: pd.DataFrame | None = None,
                   chunk: int = 200) -> dict:
    """Write what ``ExtractEdges.py`` reads, from the expression the nulls score (``adata.X``):

    - ``em.csv``: genes x cells, first column the gene symbol. Only the database's genes are written:
      NATMI reads no other gene, and its means and specificities are per gene, so nothing changes.
    - ``ann.csv``: cell, cluster. Cells are renamed ``c0, c1, ...`` so no barcode needs quoting.

    Values keep their precision: ``%.9g`` for float32 data (exact round trip), ``%.17g`` otherwise.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pairs = load_pairs() if pairs is None else pairs
    var = pd.Index(adata.var_names.astype(str))
    genes = sorted(database_genes(pairs) & set(var))
    X = as_matrix(adata.X)
    fmt = "%.9g" if getattr(adata.X, "dtype", np.float64) == np.float32 else "%.17g"
    cells = [f"c{i}" for i in range(adata.n_obs)]
    em, ann = out_dir / "em.csv", out_dir / "ann.csv"
    with open(em, "w", newline="") as f:
        f.write("gene," + ",".join(cells) + "\n")
        for a in range(0, len(genes), chunk):
            g = genes[a:a + chunk]
            block = columns(X, var.get_indexer(g))
            block = block.toarray() if hasattr(block, "toarray") else np.asarray(block)
            pd.DataFrame(block.T, index=g).to_csv(f, header=False, float_format=fmt)
    labels = adata.obs[groupby].astype(str).to_numpy()
    pd.DataFrame({"cell": cells, "cluster": labels}).to_csv(ann, index=False)
    return dict(em=em, ann=ann, out_dir=out_dir, genes=genes, n_cells=int(adata.n_obs), float_format=fmt)


def run_native(inputs: dict, out_dir: Path | None = None, *, core_num: int = 1, python: str | None = None,
               natmi: str | Path | None = None, log=None, **kw) -> dict:
    """``python ExtractEdges.py --species human --interSpecies human --idType symbol --interDB lrc2p``,
    run from the NATMI clone (it reads ``lrdbs/`` relative to its working directory). Returns the edge
    table (``Edges_lrc2p.csv``: only rows with a positive score) and how it was produced."""
    natmi = Path(natmi) if natmi is not None else natmi_dir()
    if not (natmi / "ExtractEdges.py").exists():
        raise NativeNotInstalled(f"ExtractEdges.py not found in {natmi}: clone NATMI and set NATMI_DIR "
                                 "(tools/natmi/requirements-native.txt)")
    python = python or natmi_python()
    out = (Path(out_dir) if out_dir is not None else Path(inputs["out_dir"])) / "natmi_out"
    cmd = [python, "ExtractEdges.py", "--species", "human", "--interSpecies", "human", "--idType", "symbol",
           "--emFile", str(Path(inputs["em"]).resolve()), "--annFile", str(Path(inputs["ann"]).resolve()),
           "--interDB", INTER_DB, "--coreNum", str(int(core_num)), "--out", str(out.resolve())]
    t0 = time.time()
    r = subprocess.run(cmd, cwd=natmi, capture_output=True, text=True)
    out.mkdir(parents=True, exist_ok=True)
    (out / "natmi_stdout.txt").write_text(r.stdout + "\n--- stderr ---\n" + r.stderr)
    if r.returncode != 0:
        raise RuntimeError(f"NATMI failed ({r.returncode}); see {out / 'natmi_stdout.txt'}")
    if log:
        log(f"    native NATMI: {time.time() - t0:.0f} s")
    try:
        commit = subprocess.run(["git", "-C", str(natmi), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except OSError:
        commit = None
    edges = pd.read_csv(out / f"Edges_{INTER_DB}.csv",
                        dtype={"Sending cluster": str, "Target cluster": str, "Ligand symbol": str, "Receptor symbol": str})
    return dict(edges=edges, out_dir=out, command=cmd, commit=commit, commit_expected=NATMI_COMMIT,
                seconds=round(time.time() - t0, 1))


def native_to_schema(native: dict) -> pd.DataFrame:
    """NATMI's edges as the common schema; ``score`` = edge average expression derived specificity,
    ``expr_prod`` = edge average expression weight. Only rows with a positive score (NATMI writes no others)."""
    e = native["edges"]
    lig, rec = e["Ligand symbol"].astype(str).to_numpy(), e["Receptor symbol"].astype(str).to_numpy()
    df = pd.DataFrame({"ligand": lig, "receptor": rec,
                       "source": e["Sending cluster"].astype(str).to_numpy(),
                       "target": e["Target cluster"].astype(str).to_numpy(),
                       "score": e[NATIVE_SCORE_NAME].to_numpy(float), "score_name": NATIVE_SCORE_NAME,
                       "p_native": np.nan, "ligand_complex": lig, "receptor_complex": rec, "tool": TOOL,
                       "expr_prod": e[NATIVE_WEIGHT_NAME].to_numpy(float)})
    return schema.validate(df, require_p_native=HAS_NATIVE_NULL)


def fidelity(port: pd.DataFrame, native: pd.DataFrame, rtol: float = FIDELITY_RTOL) -> dict:
    """Column 1 vs column 2 (AMENDMENT_5, T0): the rows with a positive score must be the same set, and
    ``score`` and ``expr_prod`` must agree to a relative difference of ``rtol``."""
    key = ["ligand", "receptor", "source", "target"]
    p = port.loc[port.score.to_numpy() > 0, key + ["score", "expr_prod"]]
    m = p.merge(native[key + ["score", "expr_prod"]], on=key, how="outer", suffixes=("_port", "_native"), indicator=True)
    both = m[m._merge == "both"]

    def rel(c):
        a, b = both[f"{c}_port"].to_numpy(float), both[f"{c}_native"].to_numpy(float)
        return float((np.abs(a - b) / np.maximum(np.abs(b), 1e-300)).max()) if len(both) else 0.0

    out = dict(native_rows=int(len(native)), port_positive_rows=int(len(p)),
               only_native=int((m._merge == "right_only").sum()),
               only_port=int((m._merge == "left_only").sum()), rows_compared=int(len(both)),
               max_rel_diff_score=rel("score"), max_rel_diff_expr_prod=rel("expr_prod"), rtol=rtol)
    out["identical"] = bool(out["only_native"] == 0 and out["only_port"] == 0
                            and out["max_rel_diff_score"] <= rtol and out["max_rel_diff_expr_prod"] <= rtol)
    return out
