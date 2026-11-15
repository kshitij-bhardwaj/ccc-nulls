"""CellPhoneDB (v5) adapter: dataset -> native inputs -> native run -> common schema.

Contract: tools/README.md. Common schema: ccc_nulls.schema.

The native path runs CellPhoneDB Method 1 (means only). Significance comes from this
project's Null A, computed on CellPhoneDB's own statistic (see ``cpdb_null``), not
from CellPhoneDB Method 2.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from ccc_nulls import liana_reference as _liana
from ccc_nulls import schema
from ccc_nulls.stats import bh

from . import cpdb_null

TOOL = "cellphonedb"
LIANA_METHOD = "cellphonedb"
SCORE_NAME = "lr_means"
HAS_NATIVE_NULL = True       # CellPhoneDB has a label-permutation null (Method 2)
NATIVE_MODE = "method1"      # ...but the native path here runs Method 1, so p_native is NaN
SEPARATOR = "|"              # CellPhoneDB's "<type a>|<type b>" pair-column separator
SCORE_TOLERANCE = 5e-4 + 1e-6  # CellPhoneDB rounds means to 3 decimals


def liana_reference(adata, groupby: str = "label", **kw) -> pd.DataFrame:
    """LIANA's re-implementation of CellPhoneDB (v5), in the common schema (working baseline)."""
    return _liana.run(adata, LIANA_METHOD, groupby, **kw)


def prepare_inputs(adata, groupby: str, out_dir: Path, *, layer: str | None = None,
                   counts_data: str = "hgnc_symbol", database_path=None) -> dict:
    """Write CellPhoneDB's inputs (``counts.h5ad`` + ``meta.tsv``) and return their paths.

    The expression matrix (``adata.X`` or ``adata.layers[layer]``) must be normalised and
    unscaled, with gene identifiers of type ``counts_data`` as ``var_names``.
    CellPhoneDB uses it as supplied. It is written as a sparse float32 h5ad (CellPhoneDB
    converts to float32 itself), so it is never densified here.

    With ``database_path``, only genes that appear in the CellPhoneDB database are
    written. CellPhoneDB drops every other gene when it maps the counts onto its
    database, so the results are unchanged. CellPhoneDB loads the counts densely, so
    this keeps large atlases in memory.
    """
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError("prepare_inputs needs anndata (installed with cellphonedb)") from exc
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if groupby not in adata.obs:
        raise KeyError(f"{groupby!r} is not present in adata.obs")

    cells = pd.Index(adata.obs_names).astype(str)
    genes = pd.Index(adata.var_names).astype(str)
    if not cells.is_unique or not genes.is_unique:
        raise ValueError("CellPhoneDB requires unique cell IDs and gene identifiers")
    if adata.obs[groupby].isna().any():
        raise ValueError(f"{int(adata.obs[groupby].isna().sum())} cells have no {groupby!r} label")
    labels = adata.obs[groupby].astype(str).to_numpy()
    bad = sorted({lab for lab in labels if SEPARATOR in lab})
    if bad:
        raise ValueError(f"cell-type labels must not contain {SEPARATOR!r} (CellPhoneDB's pair separator): {bad}")

    X = adata.layers[layer] if layer else adata.X
    X = sparse.csr_matrix(X, dtype=np.float32)
    if X.shape != (len(cells), len(genes)):
        raise ValueError("expression matrix dimensions do not agree with obs_names and var_names")
    if database_path is not None:
        keep = genes.isin(database_genes(database_path, counts_data))
        if not keep.any():
            raise ValueError(f"no var_names match the CellPhoneDB database's {counts_data!r} identifiers")
        X, genes = X[:, np.flatnonzero(keep)], genes[keep]

    counts_path = out_dir / "counts.h5ad"
    meta_path = out_dir / "meta.tsv"
    ad.AnnData(X=X, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes)).write_h5ad(counts_path)
    pd.DataFrame({"Cell": cells, "cell_type": labels}).to_csv(meta_path, sep="\t", index=False)
    return {"counts": counts_path, "meta": meta_path, "counts_data": counts_data}


def database_genes(database_path, counts_data: str = "hgnc_symbol") -> set:
    """Gene identifiers of type ``counts_data`` that CellPhoneDB can map onto its database."""
    try:
        from cellphonedb.utils import db_utils
    except ImportError as exc:
        raise ImportError("install cellphonedb to read its database") from exc
    genes = db_utils.get_interactions_genes_complex(str(database_path))[1]
    return set(genes[counts_data].dropna().astype(str))


def resource_pairs(database_path) -> pd.DataFrame:
    """CellPhoneDB's directed simple interactions as ``ligand, receptor`` gene symbols.

    ``database_path`` is the downloaded ``cellphonedb.zip`` or its folder. The
    ``*_input.csv`` tables that ``download_database`` writes next to it are read from
    there. Complexes and interactions where both or neither partner is a receptor are left
    out, as in :func:`native_to_schema`. The ligand is the partner that is not a receptor.
    Use it as this tool's resource for ``ccc_nulls.simulate.restrict_truth`` /
    ``plant_pairs``.
    """
    d = Path(database_path)
    d = d.parent if d.suffix == ".zip" else d
    inter = pd.read_csv(d / "interaction_input.csv")
    genes = pd.read_csv(d / "gene_input.csv").dropna(subset=["hgnc_symbol"])
    prot = pd.read_csv(d / "protein_input.csv")
    complexes = set(pd.read_csv(d / "complex_input.csv").complex_name)
    is_rec = dict(zip(prot.uniprot, _as_bool(prot.receptor)))
    simple = inter[~inter.partner_a.isin(complexes) & ~inter.partner_b.isin(complexes)]
    ra, rb = simple.partner_a.map(is_rec), simple.partner_b.map(is_rec)
    directed = simple[ra.notna() & rb.notna() & (ra != rb)]
    a_is_rec = directed.partner_a.map(is_rec).to_numpy(bool)
    ligands = np.where(a_is_rec, directed.partner_b, directed.partner_a)
    receptors = np.where(a_is_rec, directed.partner_a, directed.partner_b)
    symbols = genes.groupby("uniprot").hgnc_symbol.apply(list)
    rows = [(lig, rec) for a, b in zip(ligands, receptors)
            for lig in symbols.get(a, []) for rec in symbols.get(b, [])]
    return pd.DataFrame(rows, columns=["ligand", "receptor"]).drop_duplicates().reset_index(drop=True)


def run_native(inputs: dict, out_dir: Path, **kw) -> dict:
    """Run native CellPhoneDB Method 1 (means only, no label permutations).

    ``database_path`` is required and must point to a downloaded CellPhoneDB
    database zip.
    """
    database_path = kw.pop("database_path", None)
    if database_path is None:
        raise NotImplementedError(
            "CellPhoneDB Method 1 needs a downloaded database; pass database_path=/path/to/cellphonedb.zip"
        )
    threshold = float(kw.pop("threshold", 0.1))
    if not 0 <= threshold <= 1:
        raise ValueError("threshold must be between 0 and 1")
    if kw:
        unknown = ", ".join(sorted(kw))
        raise TypeError(f"unsupported Method 1 options: {unknown}")
    try:
        from cellphonedb.src.core.methods import cpdb_analysis_method
    except ImportError as exc:
        raise ImportError("install cellphonedb in the native-tool environment to run Method 1") from exc

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    result = cpdb_analysis_method.call(
        cpdb_file_path=str(database_path),
        meta_file_path=str(inputs["meta"]),
        counts_file_path=str(inputs["counts"]),
        counts_data=inputs.get("counts_data", "hgnc_symbol"),
        output_path=str(out_dir),
        threshold=threshold,
    )
    means = result.get("means_result", result.get("means"))
    deconvoluted = result.get("deconvoluted_result", result.get("deconvoluted"))
    if means is None:
        raise RuntimeError("CellPhoneDB Method 1 did not return a means table")
    return {"means": means, "deconvoluted": deconvoluted}


def method1_with_nulls(adata, groupby: str, out_dir: Path, *, database_path, B: int = 1000,
                       threshold: float = 0.1, seed: int = 0, pvalue: str = "ccc", alpha: float = 0.05,
                       layer: str | None = None, check_scores: bool = True) -> pd.DataFrame:
    """Run CellPhoneDB Method 1, then this project's Null A on CellPhoneDB's own statistic.

    One row per directed CellPhoneDB interaction x (source, target) cell-type pair.
    Complexes are scored the CellPhoneDB way (minimum subunit), and ``ligand`` /
    ``receptor`` name the subunit that carries that minimum in the source / target
    type. Added columns:

    ``score_obs``  CellPhoneDB mean recomputed here (``score`` is CellPhoneDB's rounded value)
    ``tested``     passed CellPhoneDB's expression gate (> ``threshold`` of cells) with score > 0
    ``p_rand``     Null A p-value (NaN where not tested)
    ``q_rand``     Benjamini-Hochberg over tested rows only (NaN elsewhere)

    Null B is not computed here yet.
    """
    out_dir = Path(out_dir)
    inputs = prepare_inputs(adata, groupby, out_dir / "inputs", layer=layer, database_path=database_path)
    native = run_native(inputs, out_dir / "method1", database_path=database_path, threshold=threshold)
    data = cpdb_null.load(inputs["counts"], inputs["meta"], database_path, inputs["counts_data"])
    null = cpdb_null.null_a(data, B=B, threshold=threshold, seed=seed, pvalue=pvalue)
    return combine(native, null, alpha=alpha, check_scores=check_scores)


def combine(native: dict, null: cpdb_null.NullAResult, *, alpha: float = 0.05,
            check_scores: bool = True) -> pd.DataFrame:
    """Join CellPhoneDB's means table with Null A results, in the common schema."""
    means, flip, pairs = _directed_means(native)
    rows = pd.Index(null.interaction_ids).get_indexer(means.id_cp_interaction)
    cols = pd.Index(null.pairs).get_indexer(pairs)
    if (rows < 0).any() or (cols < 0).any():
        raise ValueError("CellPhoneDB means and the Null A result do not cover the same "
                         "interactions / cell-type pairs; were they run on the same inputs?")
    types = pd.Index(null.cell_types)
    first, second = _split_pairs(pairs)
    gene_1 = null.gene_1[rows][:, types.get_indexer(first)]
    gene_2 = null.gene_2[rows][:, types.get_indexer(second)]
    native_score = means[pairs].to_numpy(dtype=float)
    out = _melt(means, flip, pairs, gene_1, gene_2, {
        "score": native_score,
        "score_obs": null.score[rows][:, cols],
        "tested": null.tested[rows][:, cols],
        "p_rand": null.p[rows][:, cols],
    })
    if check_scores:
        diff = (out.score - out.score_obs).abs()
        if (diff > SCORE_TOLERANCE).any():
            worst = out.loc[diff.nlargest(5).index, ["native_interaction_id", "cpdb_pair", "score", "score_obs"]]
            raise RuntimeError("recomputed scores differ from CellPhoneDB's means beyond rounding:\n"
                               f"{worst.to_string()}")
    q = np.full(len(out), np.nan)
    tested = out.tested.to_numpy(bool)
    q[tested] = bh(out.p_rand.to_numpy()[tested])
    out["q_rand"] = q
    out["sig_rand"] = out.q_rand < alpha
    return schema.validate(out)


def native_to_schema(native: dict) -> pd.DataFrame:
    """Convert Method 1's wide means table into the shared long interaction schema.

    One row per directed interaction x (source, target) pair. ``ligand`` /
    ``receptor`` are gene symbols for simple partners and the complex name for
    complexes (the means table does not say which subunit carries a complex's minimum;
    :func:`combine` fills that in). ``p_native`` is NaN because Method 1 has no p-values.
    """
    means, flip, pairs = _directed_means(native)
    gene_1 = _partner_names(means.partner_a, means.gene_a)
    gene_2 = _partner_names(means.partner_b, means.gene_b)
    n = len(pairs)
    out = _melt(means, flip, pairs, np.repeat(gene_1[:, None], n, 1), np.repeat(gene_2[:, None], n, 1),
                {"score": means[pairs].to_numpy(dtype=float)})
    finite = np.isfinite(out.score.to_numpy())
    if not finite.all():
        warnings.warn(f"dropping {int((~finite).sum())} rows with a non-numeric CellPhoneDB mean", stacklevel=2)
        out = out.loc[finite]
    return schema.validate(out)


# ------------------------------------------------------------------ helpers
def _read_table(table):
    if isinstance(table, (str, Path)):
        return pd.read_csv(table, sep=None, engine="python")
    if not isinstance(table, pd.DataFrame):
        raise TypeError("expected a DataFrame or a path to a CellPhoneDB table")
    return table


def _as_bool(values) -> np.ndarray:
    return pd.Series(values).astype(str).str.strip().str.lower().isin({"true", "1", "1.0"}).to_numpy()


def _complex_name(partner) -> np.ndarray:
    return pd.Series(partner).astype(str).str.split(":", n=1).str[-1].to_numpy()


def _partner_names(partner, gene) -> np.ndarray:
    """Gene symbol for a simple partner, complex name for a complex."""
    partner = pd.Series(partner).astype(str).to_numpy()
    is_complex = np.char.startswith(partner.astype(str), "complex:")
    return np.where(is_complex, _complex_name(partner), pd.Series(gene).astype(str).to_numpy())


def _split_pairs(pairs):
    parts = [str(p).split(SEPARATOR, 1) for p in pairs]
    return np.array([a for a, _ in parts]), np.array([b for _, b in parts])


def _directed_means(native: dict):
    """Means table restricted to ligand -> receptor interactions, the receptor flags and pair columns.

    CellPhoneDB orders partners as A in the first cell type of ``"x|y"`` and B in the
    second, not ligand first. An interaction is directed when exactly one partner is a
    receptor. ``flip`` marks rows where A is the receptor (so the ligand is B, in y).
    """
    means = _read_table(native.get("means"))
    required = {"id_cp_interaction", "partner_a", "partner_b", "gene_a", "gene_b", "receptor_a", "receptor_b"}
    missing = required - set(means.columns)
    if missing:
        raise ValueError(f"CellPhoneDB means table is missing columns {sorted(missing)}")
    pairs = [c for c in means.columns if SEPARATOR in str(c)]
    if not pairs:
        raise ValueError("CellPhoneDB means table contains no '<source>|<target>' columns")
    ra, rb = _as_bool(means.receptor_a), _as_bool(means.receptor_b)
    directed = ra != rb
    if not directed.all():
        warnings.warn(f"skipping {int((~directed).sum())} of {len(means)} CellPhoneDB interactions where both "
                      "or neither partner is a receptor (no ligand -> receptor direction)", stacklevel=3)
    if not directed.any():
        raise ValueError("no directed ligand-receptor interactions in the CellPhoneDB means table")
    return means.loc[directed].reset_index(drop=True), ra[directed], pairs


def _melt(means, flip, pairs, gene_1, gene_2, values: dict) -> pd.DataFrame:
    """Wide (interaction x pair) arrays -> long rows oriented ligand -> receptor.

    ``gene_1`` / ``gene_2`` (interactions x pairs) name the gene used for partner A in
    the pair's first type and for partner B in its second type.
    """
    n, k = len(means), len(pairs)
    ii, jj = np.repeat(np.arange(n), k), np.tile(np.arange(k), n)
    first, second = _split_pairs(pairs)
    f = np.asarray(flip, bool)[ii]
    g1, g2 = np.asarray(gene_1).ravel(), np.asarray(gene_2).ravel()
    # complex name for a complex partner; the gene itself for a simple one (the schema's
    # convention), not CellPhoneDB's "simple:<uniprot>" id
    cx_a = np.char.startswith(means.partner_a.astype(str).to_numpy().astype(str), "complex:")[ii]
    cx_b = np.char.startswith(means.partner_b.astype(str).to_numpy().astype(str), "complex:")[ii]
    n1, n2 = _complex_name(means.partner_a)[ii], _complex_name(means.partner_b)[ii]
    # CellPhoneDB has complexes named like a gene (complex "GALR1" next to gene GALR1);
    # keep its "complex:" prefix only for those, so the two stay distinct rows
    genes = set(g1[~cx_a]) | set(g2[~cx_b])
    n1 = np.where(cx_a & np.isin(n1, list(genes)), np.char.add("complex:", n1.astype(str)), n1)
    n2 = np.where(cx_b & np.isin(n2, list(genes)), np.char.add("complex:", n2.astype(str)), n2)
    c1, c2 = np.where(cx_a, n1, g1), np.where(cx_b, n2, g2)
    out = pd.DataFrame({
        "ligand": np.where(f, g2, g1), "receptor": np.where(f, g1, g2),
        "source": np.where(f, second[jj], first[jj]), "target": np.where(f, first[jj], second[jj]),
        "score_name": SCORE_NAME, "p_native": np.nan,
        "ligand_complex": np.where(f, c2, c1), "receptor_complex": np.where(f, c1, c2), "tool": TOOL,
        "native_interaction_id": means.id_cp_interaction.to_numpy()[ii],
        "cpdb_pair": np.asarray(pairs, dtype=object)[jj],
    })
    for name, wide in values.items():
        out[name] = np.asarray(wide).ravel()
    return out
