"""CellPhoneDB's own statistic, and this project's Null A (label permutation) on it.

The statistic mirrors cellphonedb 5.0.1 (``cpdb_analysis_method`` and
``cpdb_statistical_analysis_helper``):

* Expression is the counts file mapped onto CellPhoneDB multidata ids: float32,
  identifiers that map to the same id averaged, all-zero genes dropped, and only
  interactions whose partners (every subunit, for a complex) are present kept.
* Per cell type: the mean expression and the fraction of expressing cells (> 0).
* A complex takes the minimum over its subunits, for the mean and for the fraction.
* The score of partner a in cell type x with partner b in cell type y is
  ``(m_a > 0) * (m_b > 0) * (m_a + m_b) / 2``.
* Expression gate: ``frac_a > threshold and frac_b > threshold`` (strict).

Null A permutes the cell-type labels, recomputes the same score, and compares it with
the observed one. As in CellPhoneDB Method 2, permuted scores are not gated. Only
gated-in interactions with a non-zero observed score are tested. The rest get a NaN
p-value, so the caller can leave them out of multiple-testing correction.

The p-value is ``(1 + #{b : s_b >= s_obs}) / (B + 1)`` (``pvalue="ccc"``, never 0).
``pvalue="cpdb"`` gives CellPhoneDB's own ``#{b : s_b > s_obs} / B`` for
side-by-side checks. The permutations come from this module's own NumPy generator,
not from CellPhoneDB.

``load()`` is the only function that needs the ``cellphonedb`` package. It reuses
CellPhoneDB's own loading and prefiltering code, so the statistic is computed on
exactly the data CellPhoneDB sees. Everything else is plain NumPy.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


@dataclass
class CpdbData:
    """CellPhoneDB's filtered input, ready for scoring.

    ``counts`` has one row per multidata id in ``row_ids`` (simple proteins) and one
    column per cell. ``complexes`` maps a complex multidata id to the row positions of
    its subunits. ``interactions`` has columns ``multidata_1_id``, ``multidata_2_id``
    and ``id_cp_interaction``. Partner 1 sits in the first cell type of a
    ``"x|y"`` pair column and partner 2 in the second, as in CellPhoneDB's tables.
    """

    counts: sparse.csr_matrix
    row_ids: np.ndarray
    row_genes: np.ndarray
    codes: np.ndarray
    cell_types: list
    complexes: dict
    interactions: pd.DataFrame
    separator: str = "|"

    def __post_init__(self):
        self.counts = sparse.csr_matrix(self.counts, dtype=np.float64)
        self.row_ids = np.asarray(self.row_ids)
        self.row_genes = np.asarray(self.row_genes).astype(str)
        self.codes = np.asarray(self.codes, dtype=np.int64)
        self.cell_types = [str(c) for c in self.cell_types]
        if self.counts.shape != (len(self.row_ids), len(self.codes)):
            raise ValueError("counts must be rows (row_ids) x cells (codes)")
        if np.bincount(self.codes, minlength=self.n_types).min() == 0:
            raise ValueError("every cell type must contain at least one cell")
        if any(self.separator in c for c in self.cell_types):
            raise ValueError(f"cell-type names must not contain the separator {self.separator!r}")

        # Feature table: simple rows first, then complexes (the order build_clusters uses).
        self.complex_ids = list(self.complexes)
        width = max((len(r) for r in self.complexes.values()), default=1)
        pad = np.zeros((len(self.complex_ids), width), dtype=np.int64)
        for k, cid in enumerate(self.complex_ids):
            rows = np.asarray(self.complexes[cid], dtype=np.int64)
            if len(rows) == 0:
                raise ValueError(f"complex {cid!r} has no subunit rows")
            pad[k] = np.concatenate([rows, np.full(width - len(rows), rows[0])])
        self._complex_rows = pad
        position = {rid: i for i, rid in enumerate(self.row_ids)}
        position.update({cid: len(self.row_ids) + k for k, cid in enumerate(self.complex_ids)})
        missing = set(self.interactions.multidata_1_id) | set(self.interactions.multidata_2_id)
        missing -= set(position)
        if missing:
            raise ValueError(f"interaction partners missing from counts/complexes: {sorted(missing)[:5]}")
        self.i1 = self.interactions.multidata_1_id.map(position).to_numpy(np.int64)
        self.i2 = self.interactions.multidata_2_id.map(position).to_numpy(np.int64)
        self._expressed = (self.counts > 0).astype(np.float64).tocsr()

    @property
    def n_types(self) -> int:
        return len(self.cell_types)

    @property
    def pairs(self) -> list:
        """Pair columns in CellPhoneDB's order: ``x|y`` for x, then y, over sorted types."""
        s = self.separator
        return [f"{a}{s}{b}" for a in self.cell_types for b in self.cell_types]


# ------------------------------------------------------------------ statistic
def _type_average(M: sparse.csr_matrix, codes: np.ndarray, n_types: int) -> np.ndarray:
    """rows x types mean of M's columns within each cell type."""
    onehot = np.zeros((len(codes), n_types))
    onehot[np.arange(len(codes)), codes] = 1.0
    return np.asarray(M @ onehot) / onehot.sum(0)


def feature_stats(data: CpdbData, codes=None, *, percents: bool = True):
    """Per-cell-type mean (and fraction expressing) for every simple protein and complex.

    Returns ``(means, fracs, rep)``, each features x types. ``rep`` is the row (simple
    protein) that carries a complex's minimum mean in that cell type, which is the
    subunit CellPhoneDB's score effectively uses.
    """
    codes = data.codes if codes is None else codes
    m = _type_average(data.counts, codes, data.n_types)
    sub = m[data._complex_rows]                              # complexes x width x types
    arg = sub.argmin(1)
    means = np.vstack([m, np.take_along_axis(sub, arg[:, None, :], 1)[:, 0, :]])
    rep = np.vstack([np.repeat(np.arange(len(m))[:, None], data.n_types, 1),
                     np.take_along_axis(data._complex_rows[:, :, None].repeat(data.n_types, 2),
                                        arg[:, None, :], 1)[:, 0, :]])
    if not percents:
        return means, None, rep
    f = _type_average(data._expressed, codes, data.n_types)
    fracs = np.vstack([f, f[data._complex_rows].min(1)])
    return means, fracs, rep


def pair_scores(means: np.ndarray, i1, i2) -> np.ndarray:
    """interactions x (types * types) CellPhoneDB mean score; column ``a * K + b`` is ``a|b``."""
    x, y = means[i1], means[i2]
    s = ((x > 0)[:, :, None] & (y > 0)[:, None, :]) * (x[:, :, None] + y[:, None, :]) / 2
    return s.reshape(len(x), -1)


def pair_gate(fracs: np.ndarray, i1, i2, threshold: float) -> np.ndarray:
    """interactions x (types * types): both partners expressed in > threshold of their cells."""
    x, y = fracs[i1] > threshold, fracs[i2] > threshold
    return (x[:, :, None] & y[:, None, :]).reshape(len(x), -1)


# ------------------------------------------------------------------ Null A
@dataclass
class NullAResult:
    """Wide results, rows = ``interaction_ids`` and columns = ``pairs``."""

    interaction_ids: np.ndarray
    pairs: list
    cell_types: list
    score: np.ndarray          # observed CellPhoneDB mean (unrounded)
    tested: np.ndarray         # passed the expression gate and score > 0
    p: np.ndarray              # Null A p-value, NaN where not tested
    gene_1: np.ndarray         # interactions x types: gene carrying partner 1 in each type
    gene_2: np.ndarray         # interactions x types: gene carrying partner 2 in each type
    B: int
    pvalue: str
    threshold: float
    seed: int
    meta: dict = field(default_factory=dict)

    def wide(self, name: str) -> pd.DataFrame:
        """One result as a CellPhoneDB-shaped table (``score``, ``tested`` or ``p``)."""
        return pd.DataFrame(getattr(self, name), index=pd.Index(self.interaction_ids, name="id_cp_interaction"),
                            columns=self.pairs)


def null_a(data: CpdbData, *, B: int = 1000, threshold: float = 0.1, seed: int = 0,
           pvalue: str = "ccc") -> NullAResult:
    """Uniform label-permutation null on CellPhoneDB's own statistic."""
    if pvalue not in ("ccc", "cpdb"):
        raise ValueError("pvalue must be 'ccc' ((1 + #>=) / (B + 1)) or 'cpdb' (#> / B)")
    if B < 1:
        raise ValueError("B must be at least 1")
    means, fracs, rep = feature_stats(data)
    s_obs = pair_scores(means, data.i1, data.i2)
    tested = pair_gate(fracs, data.i1, data.i2, threshold) & (s_obs > 0)

    rng = np.random.default_rng(seed)
    exceed = np.zeros(s_obs.shape, dtype=np.int64)
    for _ in range(B):
        m_perm, _, _ = feature_stats(data, rng.permutation(data.codes), percents=False)
        s_perm = pair_scores(m_perm, data.i1, data.i2)
        exceed += (s_perm > s_obs) if pvalue == "cpdb" else (s_perm >= s_obs)
    p = exceed / B if pvalue == "cpdb" else (1 + exceed) / (B + 1)
    p = np.where(tested, p, np.nan)

    genes = data.row_genes
    return NullAResult(
        interaction_ids=data.interactions.id_cp_interaction.to_numpy(), pairs=data.pairs,
        cell_types=list(data.cell_types), score=s_obs, tested=tested, p=p, gene_1=genes[rep[data.i1]], gene_2=genes[rep[data.i2]],
        B=B, pvalue=pvalue, threshold=threshold, seed=seed,
    )


# ------------------------------------------------------------------ loading
def load(counts, meta_path, database_path, counts_data: str = "hgnc_symbol") -> CpdbData:
    """Load and prefilter input exactly as ``cpdb_analysis_method.call`` (cellphonedb 5.0.x) does.

    ``counts`` is a counts file path (``.h5ad``, ``.txt``/``.tsv``, ...) or an AnnData
    object, as accepted by CellPhoneDB itself. This reuses CellPhoneDB's private
    helpers, so it is pinned to the 5.0 series. Check it again when upgrading.
    """
    try:
        from cellphonedb.src.core.methods import cpdb_statistical_analysis_helper as helper
        from cellphonedb.src.core.models.complex import complex_helper
        from cellphonedb.utils import db_utils, file_utils
    except ImportError as exc:
        raise ImportError("install cellphonedb (5.0.x) to load CellPhoneDB data") from exc

    interactions, genes, complex_compositions, complexes, gene_synonym2gene_name, _ = \
        db_utils.get_interactions_genes_complex(str(database_path))
    counts = counts if not isinstance(counts, (str, Path)) else str(counts)
    counts, meta, _, _, _ = file_utils.get_user_files(
        counts=counts, meta_fp=str(meta_path), gene_synonym2gene_name=gene_synonym2gene_name,
        counts_data=counts_data)
    reduced = interactions[["multidata_1_id", "multidata_2_id"]].drop_duplicates()
    counts, relations = helper.add_multidata_and_means_to_counts(counts, genes, counts_data)
    if counts.empty:
        raise ValueError("no counts rows map onto the CellPhoneDB database (are these human gene symbols?)")
    filtered, counts_filtered, composition = helper.prefilters(reduced, counts, complexes, complex_compositions)
    if filtered.empty:
        raise ValueError("no CellPhoneDB interaction has all of its partners expressed in this input")

    meta = meta.loc[counts_filtered.columns]
    labels = pd.Categorical(meta["cell_type"].astype(str))
    complex_rows = complex_helper.map_complex_to_protein_row_ids(composition, counts_filtered)
    symbol = relations.drop_duplicates("id_multidata").set_index("id_multidata")[counts_data]
    filtered = filtered.copy()
    filtered["id_cp_interaction"] = interactions.loc[filtered.index, "id_cp_interaction"].to_numpy()
    return CpdbData(
        counts=sparse.csr_matrix(counts_filtered.to_numpy(np.float32)),
        row_ids=counts_filtered.index.to_numpy(),
        row_genes=symbol.reindex(counts_filtered.index).astype(str).to_numpy(),
        codes=labels.codes, cell_types=list(labels.categories),
        complexes={cid: list(rows) for cid, rows in complex_rows.items()},
        interactions=filtered.reset_index(drop=True),
    )
