"""The tools, behind one interface: native run (column 1), our Null A (column 2), Null B (column 3).

A tool here is a thin wrapper over its own adapter (``tools/<tool>/adapter.py``); nothing in the
tool folders is changed. ``Tool.run(adata, groupby, workdir)`` returns one row per
(L-R, source, target) with:

    score             the tool's own statistic (CellChat prob, CellPhoneDB mean)
    p_native          column 1: the tool's own p-value, from the tool itself
    p_port_native     column 2: our Null A, same permutations and p formula as the tool
    call_native, call_port_native   the tool's own call rule on each
    p_null_a          our Null A with B permutations and the repository's p convention
    p_null_b, p_L_null_b, p_R_null_b, limiting_partner   Null B, same permutations as p_null_a
    family            the rows BH corrects over (the family frozen for each tool), same for both nulls
    ligand_complex, receptor_complex   subunit symbols joined by "_", for joins with LIANA+

and ``meta`` with the fidelity check (column 1 vs column 2) and timings. Null A and Null B use the
tool's own permutation sequence: CellChat's (R's ``set.seed(1); sample.int``), CellPhoneDB's
(``np.random.seed(0)`` then one ``shuffle`` per iteration).
"""
from __future__ import annotations

import contextlib
import io
import time
from pathlib import Path

import numpy as np
import pandas as pd

from ccc_nulls import composite


def _secs(t0):
    return round(time.time() - t0, 1)


def complex_key(names) -> np.ndarray:
    """Subunits sorted and joined by "_", so complexes match across tools and LIANA+."""
    return np.array(["_".join(sorted(str(x).split("_"))) for x in names], dtype=object)


class CellChat:
    name = "cellchat"
    liana_method = "cellchat"
    liana_resource = "cellchatdb"
    native_rule = "p < 0.05 and prob > 0 (subsetCommunication, thresh = 0.05)"

    @staticmethod
    def tool_rule(p, df) -> np.ndarray:
        """CellChat's call rule applied to any p-value column."""
        return (np.asarray(p, float) < 0.05) & (df.score.to_numpy(float) > 0)
    native_B = 100          # computeCommunProb's nboot
    native_seed = 1         # computeCommunProb's seed.use

    def __init__(self, B: int = 1000, r_cache: Path | None = None, run_r: bool = True, species: str = "human"):
        self.B, self.r_cache, self.do_r, self.species = B, r_cache, run_r, species
        if species != "human":
            self.liana_resource = f"{species}consensus"     # LIANA+ ships CellChatDB for human only

    def run(self, adata, groupby: str, workdir: Path, log=print) -> tuple[pd.DataFrame, dict]:
        from tools.cellchat import adapter, native_r, nulls
        from tools.cellchat import cellchat_py as cc

        meta, t0 = {}, time.time()
        inputs = adapter.prepare_inputs(adata, groupby, species=self.species)
        sig = None if self.species == "human" else cc.extract_gene(cc.subset_db(cc.load_db(self.species)))
        port = adapter.run_native(inputs, nboot=2, seed=self.native_seed)      # rows and prob; p from null_a below
        df = adapter.native_to_schema(port)
        meta["seconds_prepare"] = _secs(t0)

        # column 2 at the tool's own setting: CellChat's permutations and p formula (#{prob_b > prob} / nboot)
        t0 = time.time()
        a_nat = nulls.null_a(adata, df, groupby, B=self.native_B, seed=self.native_seed, rng="cellchat",
                             genes_signaling=sig)
        df["p_port_native"] = a_nat.p_null_a_cellchat.to_numpy()
        meta["seconds_port_native"] = _secs(t0)

        # column 1: CellChat itself, in R
        if self.do_r:
            t0 = time.time()
            r = native_r.run_r(adata, groupby, self.r_cache, nboot=self.native_B, seed=self.native_seed,
                               species=self.species, log=log)
            meta["seconds_native"] = _secs(t0)
            meta["native_session"] = r["session"]
            port_full = dict(prob=None, LRsig=port["LRsig"], groups=port["groups"], features_sig=inputs["features_sig"])
            port_full["prob"] = a_nat.score_null_a.to_numpy().reshape(len(port["LRsig"]), len(port["groups"]), -1)
            port_full["pval"] = df.p_port_native.to_numpy().reshape(port_full["prob"].shape)
            meta["fidelity"] = native_r.fidelity(r, port_full)
            df["p_native"] = r["pval"].ravel() if meta["fidelity"]["lr_identical"] else np.nan
            df["score_native"] = r["prob"].ravel() if meta["fidelity"]["lr_identical"] else np.nan
        else:
            df["p_native"] = np.nan
        for col, p in (("call_native", df.p_native), ("call_port_native", df.p_port_native)):
            df[col] = (p.to_numpy() < 0.05) & (df.score.to_numpy() > 0)

        # Null A and Null B at B permutations, CellChat's permutation sequence extended to B
        t0 = time.time()
        a = nulls.null_a(adata, df, groupby, B=self.B, seed=self.native_seed, rng="cellchat", genes_signaling=sig)
        meta["seconds_null_a"] = _secs(t0)
        t0 = time.time()
        b = adapter.null_b(adata, df, groupby, B=self.B, seed=self.native_seed, rng="cellchat", genes_signaling=sig)
        meta["seconds_null_b"] = _secs(t0)
        if not np.array_equal(b.score_null_b.to_numpy(), a.score_null_a.to_numpy()):
            raise RuntimeError("Null B scored something other than CellChat's prob")
        df["p_null_a"] = a.p_null_a.to_numpy()
        for c in ("p_null_b", "p_L_null_b", "p_R_null_b"):
            df[c] = b[c].to_numpy()
        df["limiting_partner"] = np.where(b.limiting_side_null_b.to_numpy() == "L", "ligand", "receptor")
        df["family"] = adapter.testable(a, b)
        df["ligand_complex"], df["receptor_complex"] = complex_key(df.ligand_complex), complex_key(df.receptor_complex)
        meta.update(B=self.B, seed=self.native_seed, rng="cellchat (R's generator)", species=self.species, n_rows=len(df),
                    n_family=int(df.family.sum()), n_lr=int(len(port["LRsig"])), n_types=len(port["groups"]),
                    family_rule="composite.testable: some labelling (observed, Null A or either Null B side) gives prob > 0")
        return df, meta


def cellphonedb_side_score(data) -> composite.SideScore:
    """CellPhoneDB 5.0.1's own mean score, computed its way (``build_clusters``, ``mean_analysis``):
    per-type means by ``numpy_groupies.aggregate`` on the float32 counts (float32 out), a complex's
    minimum over its subunits, then ``(x > 0) * (y > 0) * (x + y) / 2``. The float64 statistic of
    ``cpdb_null.feature_stats`` agrees except at float32 near-ties, where it moves a p-value by 1/B.
    Partner a is the first type of "x|y", partner b the second (as ``null_b_cellphonedb.side_score``)."""
    import numpy_groupies as npg

    from tools.cellphonedb import cpdb_null

    counts = np.asarray(data.counts.todense(), dtype=np.float32)        # proteins x cells, as CellPhoneDB holds them

    def summarize(codes, K):
        m = npg.aggregate(np.asarray(codes), counts, func="mean", axis=1, size=K)
        return np.vstack([m, m[data._complex_rows].min(1)])

    def score(m_a, m_b):
        return cpdb_null.pair_scores(np.vstack([m_a, m_b]), data.i1, data.i2 + len(m_a))

    return composite.SideScore(summarize=summarize, score=score, name="cellphonedb 5.0.1 mean (float32, its own code)")


class CellPhoneDB:
    name = "cellphonedb"
    liana_method = "cellphonedb"
    liana_resource = "cellphonedb"
    native_rule = "p <= 0.05 (significant_means, pvalue = 0.05)"

    @staticmethod
    def tool_rule(p, df) -> np.ndarray:
        """CellPhoneDB's call rule applied to any p-value column (rows outside its gate have p = 1)."""
        return (np.asarray(p, float) <= 0.05) & df.family.to_numpy(bool)
    native_seed = 0         # debug_seed; threads = 1 so the permutation sequence is fixed
    threshold = 0.1

    def __init__(self, database_path, B: int = 1000, run_native: bool = True, species: str = "human"):
        self.db, self.B, self.do_native, self.species = str(database_path), B, run_native, species
        if species != "human":
            self.liana_resource = f"{species}consensus"
        self.native_B = B   # CellPhoneDB's default iterations = 1000: one run serves columns 1-3

    def permutations(self, n_cells: int, B: int) -> np.ndarray:
        """CellPhoneDB's own sequence: ``np.random.seed(seed)`` then, per iteration, a Fisher-Yates
        shuffle of the cell-type column (legacy generator). Shuffling 0..n-1 the same way gives the
        cell indices, cells x B."""
        rs = np.random.RandomState(self.native_seed)
        out = np.empty((n_cells, B), dtype=np.int32)
        for b in range(B):
            idx = np.arange(n_cells, dtype=np.int32)
            rs.shuffle(idx)
            out[:, b] = idx
        return out

    def run(self, adata, groupby: str, workdir: Path, log=print) -> tuple[pd.DataFrame, dict]:
        from cellphonedb.src.core.methods import cpdb_statistical_analysis_method as method2

        from tools.cellphonedb import adapter, cpdb_null

        meta, t0 = {}, time.time()
        workdir = Path(workdir)
        if self.species != "human":         # CellPhoneDB is human only: genes renamed to 1:1 human orthologs
            from .data import to_human
            adata, meta["orthologs"] = to_human(adata, self.species)
        inputs = adapter.prepare_inputs(adata, groupby, workdir / "inputs", database_path=self.db)
        data = cpdb_null.load(inputs["counts"], inputs["meta"], self.db, inputs["counts_data"])
        K = data.n_types
        meta["seconds_prepare"] = _secs(t0)

        # column 1: CellPhoneDB's statistical analysis (Method 2) itself; without it, Method 1 (means only)
        t0 = time.time()
        ids = data.interactions.id_cp_interaction.to_numpy()
        with contextlib.redirect_stdout(io.StringIO()):
            if self.do_native:
                native = method2.call(cpdb_file_path=self.db, meta_file_path=str(inputs["meta"]),
                                      counts_file_path=str(inputs["counts"]), counts_data=inputs["counts_data"],
                                      output_path=str(workdir / "method2"), iterations=self.native_B,
                                      threshold=self.threshold, threads=1, debug_seed=self.native_seed)
                P_native = native["pvalues"].set_index("id_cp_interaction").loc[ids, data.pairs].to_numpy(float)
            else:
                native = adapter.run_native(inputs, workdir / "method1", database_path=self.db, threshold=self.threshold)
                P_native = None
        meta["seconds_native"] = _secs(t0)

        # columns 2 and the nulls, on CellPhoneDB's permutation sequence
        t0 = time.time()
        perms = self.permutations(len(data.codes), self.B)
        _, fracs, rep = cpdb_null.feature_stats(data)
        stat = cellphonedb_side_score(data)

        def pair_score(codes):
            m = stat.summarize(codes, K)
            return cpdb_null.pair_scores(m, data.i1, data.i2)

        s_obs = pair_score(data.codes)
        tested = cpdb_null.pair_gate(fracs, data.i1, data.i2, self.threshold) & (s_obs > 0)
        k_gt = np.zeros(s_obs.shape, np.int32)
        k_ge = np.zeros(s_obs.shape, np.int32)
        for b in range(self.B):
            s = pair_score(data.codes[perms[:, b]])
            k_gt += s > s_obs
            k_ge += s >= s_obs
        p_port_native = np.where(tested, k_gt / self.B, 1.0)          # CellPhoneDB's #{s_b > s} / B
        p_null_a = np.where(tested, (1 + k_ge) / (self.B + 1), np.nan)  # repository convention, as cpdb_null.null_a
        meta["seconds_null_a"] = _secs(t0)
        if P_native is None:
            P_native = np.full(p_port_native.shape, np.nan)
        else:
            same = p_port_native == P_native
            meta["fidelity"] = dict(n_rows=int(same.size), pval_identical=int(same.sum()),
                                    pval_max_abs_diff=float(np.abs(p_port_native - P_native).max()),
                                    calls_native=int((P_native <= 0.05).sum()), calls_python=int((p_port_native <= 0.05).sum()),
                                    calls_identical=bool(((P_native <= 0.05) == (p_port_native <= 0.05)).all()))

        t0 = time.time()
        nb = composite.composite(stat, data.codes, K, permutations=perms, ties="count")
        meta["seconds_null_b"] = _secs(t0)
        if not np.array_equal(nb["score"], s_obs):
            raise RuntimeError("Null B scored something other than CellPhoneDB's mean")

        # long format through the adapter's own join (orientation: which partner is the ligand)
        null = cpdb_null.NullAResult(interaction_ids=ids, pairs=data.pairs, cell_types=list(data.cell_types),
                                     score=s_obs, tested=tested, p=p_null_a, gene_1=data.row_genes[rep[data.i1]],
                                     gene_2=data.row_genes[rep[data.i2]], B=self.B, pvalue="ccc",
                                     threshold=self.threshold, seed=self.native_seed)
        df = adapter.combine({"means": native["means"]}, null, check_scores=False)
        diff = (df.score - df.score_obs).abs().to_numpy()
        if (diff > 5e-4 + 2e-5 * df.score_obs.abs().to_numpy() + 1e-6).any():      # 3-decimal output, float32
            raise RuntimeError(f"recomputed scores differ from CellPhoneDB's means: max {diff.max():.2e}")
        r = pd.Index(ids).get_indexer(df.native_interaction_id)
        c = pd.Index(data.pairs).get_indexer(df.cpdb_pair)
        df["score"] = df.score_obs                      # unrounded, as both nulls test it
        df["p_native"], df["p_port_native"] = P_native[r, c], p_port_native[r, c]
        df["call_native"], df["call_port_native"] = df.p_native <= 0.05, df.p_port_native <= 0.05
        df["p_null_a"] = df.p_rand.fillna(1.0)
        for k, v in composite.to_columns(nb, (r, c)).items():
            if k in ("p_null_b", "p_L_null_b", "p_R_null_b", "limiting_side_null_b", "n_pos_L_null_b", "n_pos_R_null_b"):
                df[k] = v
        a_is_ligand = df.source.to_numpy().astype(str) == np.array([p.split("|", 1)[0] for p in df.cpdb_pair])
        df["limiting_partner"] = np.where((df.limiting_side_null_b.to_numpy() == "L") == a_is_ligand, "ligand", "receptor")
        df["family"] = df.tested.to_numpy(bool)
        df["ligand_complex"], df["receptor_complex"] = complex_key(df.ligand_complex), complex_key(df.receptor_complex)
        meta.update(B=self.B, seed=self.native_seed, rng="cellphonedb (np.random.seed, shuffle)", species=self.species,
                    n_rows=len(df),
                    n_family=int(df.family.sum()), n_lr=int(len(ids)), n_types=K,
                    family_rule="CellPhoneDB's gate: both partners in > 10% of their type's cells, and score > 0",
                    score_max_abs_diff=float(diff.max()))
        return df, meta
