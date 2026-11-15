"""Null B (the shared composite null, ccc_nulls.composite) on CellPhoneDB, through the CellPhoneDB adapter.

    python -m benchmarks.null_b_cellphonedb --db <cellphonedb.zip> --part calibration --tissue mammary
    python -m benchmarks.null_b_cellphonedb --db <cellphonedb.zip> --part sensitivity --tissue mammary
    python -m benchmarks.null_b_cellphonedb --db <cellphonedb.zip> --part tnbc

Needs ``pip install cellphonedb==5.0.1`` and the database (``db_utils.download_database(<dir>, "v5.0.0")``).
Results and how they meet the checklist: benchmarks/results/null_b_cellphonedb/README.md.

Nothing in tools/cellphonedb is changed. Every data set goes through that folder's adapter as it
is: ``prepare_inputs`` -> ``run_native`` (CellPhoneDB 5.0.1 Method 1) -> ``cpdb_null.load`` ->
``cpdb_null.null_a`` -> ``combine`` (common schema; Null A p, q over the tested rows; the check
that every recomputed score equals CellPhoneDB's own mean).

Null B is the shared composite null only (``composite.composite``: p = max(p_a, p_b), each side's
labels shuffled alone). The one CellPhoneDB-specific piece is ``side_score``, built from the
adapter's own functions: per-type means with a complex at its minimum subunit
(``cpdb_null.feature_stats``) and the score ``(m_a > 0)(m_b > 0)(m_a + m_b) / 2``
(``cpdb_null.pair_scores``), partner a read from one labelling and partner b from another. Same B,
seed, permutations and tie rule (s_b >= s) as Null A: only the null differs.

Families and q-values (BH; calls at q <= alpha):
- null_a, null_b: BH over the adapter's ``tested`` rows (CellPhoneDB's expression gate, both
  partners in > 10% of their type's cells, and score > 0). Same rows for both nulls. Primary.
- null_b|testable: BH over ``composite.testable`` rows (some labelling gives score > 0), calls kept
  where the gate passes. Sensitivity: the gate is chosen on the observed labels.
- null_a|by, null_b|by: Benjamini-Yekutieli over the tested rows (valid under any dependence).

Parts:
- calibration (Tabula Sapiens, all cells, CellPhoneDB genes): real labels; receptor-side and
  ligand-side genes moved to random cells by one permutation (every pair whose receptor / ligand
  partner is entirely shuffled is null for "s signals to t"); fully shuffled labels (every call is
  false). FDR = mean over replicates of V / R (a lower bound); fully shuffled: P(any call) with
  Clopper-Pearson. Replicate 0 of each one-sided scenario also gives type-I error and p quantiles.
- sensitivity: one receptor-shuffle data set at B = 500, 1,000, 5,000.
- tnbc (held out, run once): calls, overlap, A-only by limiting partner and expression decile, Open
  Problems cytokine truth with bootstrap CIs, fake gene pairs, patient-half reproducibility.

Writes benchmarks/results/null_b_cellphonedb/<part>_<tissue>.json.
"""
from __future__ import annotations

import argparse
import json
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy import stats as sstats

from ccc_nulls import composite, datasets
from ccc_nulls.nulls import gene_strata
from ccc_nulls.stats import bh
from tools.cellphonedb import adapter, cpdb_null

HERE = Path(__file__).resolve().parent
OUT = HERE / "results" / "null_b_cellphonedb"
ALPHA, ALPHAS, THRESHOLD = 0.05, (0.01, 0.05, 0.10), 0.1
SEEDS = {"receptor": 1000, "ligand": 2000, "labels": 3000}
QS = np.linspace(0, 1, 101)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def by(p):
    """Benjamini-Yekutieli q-values."""
    return np.minimum(1.0, bh(p) * np.sum(1.0 / np.arange(1, len(p) + 1))) if len(p) else p


def over(p, rows, fn=bh):
    q = np.ones(len(p))
    if rows.any():
        q[rows] = fn(p[rows])
    return q


# ------------------------------------------------------------------ the CellPhoneDB side score
def side_score(data: cpdb_null.CpdbData) -> composite.SideScore:
    """CellPhoneDB's score split into partner a (first type of "x|y") and partner b (second type)."""
    def summarize(codes, K):
        return cpdb_null.feature_stats(data, codes, percents=False)[0]

    def score(m_a, m_b):
        return cpdb_null.pair_scores(np.vstack([m_a, m_b]), data.i1, data.i2 + len(m_a))

    return composite.SideScore(summarize=summarize, score=score, name="cellphonedb 5.0.1 mean")


def nulls_on(data, B, seed):
    """Null A (the adapter's) and Null B (shared) on one CpdbData, same permutations."""
    t0 = time.time()
    na = cpdb_null.null_a(data, B=B, threshold=THRESHOLD, seed=seed)
    t1 = time.time()
    nb = composite.composite(side_score(data), data.codes, data.n_types, B=B, seed=seed)
    t2 = time.time()
    assert np.array_equal(nb["score"], na.score), "Null B scored something other than CellPhoneDB's mean"
    return na, nb, dict(seconds_null_a=round(t1 - t0, 1), seconds_null_b=round(t2 - t1, 1))


def procedures(p_a, p_b, tested, fam):
    return {"null_a": over(p_a, tested), "null_b": over(p_b, tested),
            "null_b|testable": np.where(tested, over(p_b, fam), 1.0),
            "null_a|by": over(p_a, tested, by), "null_b|by": over(p_b, tested, by)}


def run(adata, groupby, db, B, seed, workdir):
    """The adapter end to end (Null A), then Null B; one row per directed interaction x (source, target)."""
    inputs = adapter.prepare_inputs(adata, groupby, workdir / "inputs", database_path=db)
    native = adapter.run_native(inputs, workdir / "method1", database_path=db, threshold=THRESHOLD)
    data = cpdb_null.load(inputs["counts"], inputs["meta"], db, inputs["counts_data"])
    na, nb, secs = nulls_on(data, B, seed)
    # combine's own check allows 5.01e-4 (3-decimal rounding); CellPhoneDB averages in float32, so
    # values on a rounding boundary differ by a little more. Same check, float32-aware tolerance.
    res = adapter.combine(native, na, alpha=ALPHA, check_scores=False)
    diff = (res.score - res.score_obs).abs().to_numpy()
    if (diff > 5e-4 + 2e-5 * res.score_obs.abs().to_numpy() + 1e-6).any():
        raise RuntimeError(f"recomputed scores differ from CellPhoneDB's means: max {diff.max():.2e}")
    secs.update(score_max_abs_diff=float(diff.max()), n_over_adapter_tolerance=int((diff > adapter.SCORE_TOLERANCE).sum()))
    r = pd.Index(na.interaction_ids).get_indexer(res.native_interaction_id)
    c = pd.Index(na.pairs).get_indexer(res.cpdb_pair)
    for k, v in composite.to_columns(nb, (r, c)).items():
        res[k] = v
    # side "L" of the engine is partner a (first type of the pair); map it to ligand / receptor
    a_is_ligand = res.source.to_numpy().astype(str) == np.array([p.split("|", 1)[0] for p in res.cpdb_pair])
    lim = res.limiting_side_null_b.to_numpy()
    res["limiting_partner_null_b"] = np.where((lim == "L") == a_is_ligand, "ligand", "receptor")
    tested = res.tested.to_numpy(bool)
    fam = composite.testable(res.score_null_b, res.n_pos_L_null_b, res.n_pos_R_null_b)
    q = procedures(res.p_rand.fillna(1.0).to_numpy(), res.p_null_b.to_numpy(), tested, fam)
    meta = dict(n_rows=int(len(res)), n_tested=int(tested.sum()), n_testable=int(fam.sum()),
                n_interactions=int(res.native_interaction_id.nunique()), n_types=int(data.n_types), **secs)
    return res, q, meta, data


# ------------------------------------------------------------------ data
def load_tissue(tissue: str, db):
    """Tabula Sapiens, all cells (``datasets.load_tabula_sapiens``), CellPhoneDB genes only (CellPhoneDB
    drops every other gene; library sizes still use all genes), "|" in labels replaced (CellPhoneDB's
    pair separator)."""
    adata = datasets.load_tabula_sapiens(f"tabula_sapiens_{tissue}", genes=adapter.database_genes(db))
    adata.obs["label"] = adata.obs.label.astype(str).str.replace("|", "/", regex=False).astype("category")
    return adata


def side_genes(db) -> tuple[set, set, dict]:
    """CellPhoneDB genes only ever in a ligand partner / only ever in a receptor partner (directed
    interactions: exactly one partner is a receptor), and the genes of every complex."""
    d = Path(db).parent if str(db).endswith(".zip") else Path(db)
    inter = pd.read_csv(d / "interaction_input.csv")
    prot = pd.read_csv(d / "protein_input.csv")
    genes = pd.read_csv(d / "gene_input.csv").dropna(subset=["hgnc_symbol"])
    cpx = pd.read_csv(d / "complex_input.csv")
    sym = genes.groupby("uniprot").hgnc_symbol.apply(set).to_dict()
    is_rec = dict(zip(prot.uniprot, adapter._as_bool(prot.receptor)))
    ucols = [c for c in cpx.columns if c.startswith("uniprot_")]
    members, cpx_rec = {}, {}
    for _, row in cpx.iterrows():
        us = [u for u in row[ucols] if isinstance(u, str)]
        members[row.complex_name] = set().union(*(sym.get(u, set()) for u in us))
    cpx_rec = dict(zip(cpx.complex_name, adapter._as_bool(cpx.receptor)))

    def partner(x):
        return (members[x], cpx_rec[x]) if x in members else (sym.get(x, set()), is_rec.get(x))

    lig, rec = set(), set()
    for _, row in inter.iterrows():
        (ga, ra), (gb, rb) = partner(row.partner_a), partner(row.partner_b)
        if ra is None or rb is None or bool(ra) == bool(rb):
            continue
        lig |= gb if ra else ga
        rec |= ga if ra else gb
    return lig - rec, rec - lig, members


def known_null(res, members, shuffled, side):
    """Rows whose ligand (side = "ligand_complex") or receptor partner has every gene shuffled."""
    def genes(x):
        x = str(x).removeprefix("complex:")
        return members.get(x, {x})
    return res[side].map(lambda x: genes(x) <= shuffled).to_numpy(bool)


def shuffle_genes(adata, genes: set, rng):
    """Copy of ``adata`` with the columns of ``genes`` moved to other cells by one random permutation."""
    m = adata.var_names.isin(list(genes)).astype(float)
    X = sparse.csr_matrix(adata.X)
    out = adata.copy()
    out.X = sparse.csr_matrix(X @ sparse.diags(1 - m) + X[rng.permutation(adata.n_obs)] @ sparse.diags(m))
    return out


def shuffle_labels(adata, rng):
    out = adata.copy()
    out.obs["label"] = pd.Categorical(rng.permutation(adata.obs["label"].astype(str).to_numpy()))
    return out


def clopper_pearson(k, n, a=0.05):
    lo = sstats.beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0
    hi = sstats.beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return [float(lo), float(hi)]


def mean_ci(x):
    x = np.asarray(x, float)
    if len(x) < 2:
        return dict(mean=float(x.mean()), ci95=None)
    h = sstats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return dict(mean=float(x.mean()), ci95=[float(x.mean() - h), float(x.mean() + h)])


def describe(adata, res, q):
    """Calls, overlap and the A-only rows: which partner limits Null B, and expression decile."""
    a, b = q["null_a"] <= ALPHA, q["null_b"] <= ALPHA
    lim = res.limiting_partner_null_b.to_numpy()
    out = dict(calls={m: {str(al): int((v <= al).sum()) for al in ALPHAS} for m, v in q.items()},
               overlap=dict(both=int((a & b).sum()), null_a_only=int((a & ~b).sum()), null_b_only=int((~a & b).sum())))
    ao = a & ~b
    out["null_a_only_limiting_partner"] = {k: int((ao & (lim == k)).sum()) for k in ("ligand", "receptor")}
    mu_all = pd.Series(np.asarray(adata.X.mean(0)).ravel(), index=adata.var_names)
    for side in ("ligand", "receptor"):
        mu = res[side].map(mu_all).to_numpy(float)
        tested = res.tested.to_numpy(bool) & np.isfinite(mu)
        dec = np.full(len(res), -1)
        dec[tested] = pd.qcut(pd.Series(mu[tested]).rank(method="first"), 10, labels=False).to_numpy()
        out[f"by_{side}_decile"] = [dict(decile=d, rows=int((dec == d).sum()), null_a=int((a & (dec == d)).sum()),
                                         null_b=int((b & (dec == d)).sum()), a_only=int((ao & (dec == d)).sum()))
                                    for d in range(10)]
    return out


# ------------------------------------------------------------------ parts
def part_calibration(args):
    adata = load_tissue(args.tissue, args.db)
    lig_only, rec_only, per = side_genes(args.db)
    log(f"{args.tissue}: {adata.n_obs:,} cells, {adata.n_vars} genes, {adata.obs.label.nunique()} types; "
        f"ligand-only {len(lig_only)}, receptor-only {len(rec_only)} genes")
    out = dict(tissue=args.tissue, n_cells=int(adata.n_obs), n_types=int(adata.obs.label.nunique()), B=args.B,
               alpha=ALPHA, threshold=THRESHOLD, reps=args.reps, label_reps=args.label_reps, scenarios={})

    if not args.skip_real:
        with tempfile.TemporaryDirectory() as tmp:
            res, q, meta, _ = run(adata, "label", args.db, args.B, seed=0, workdir=Path(tmp))
        out["real"] = dict(meta, **describe(adata, res, q))
        log("real:", out["real"]["calls"], out["real"]["overlap"], out["real"]["null_a_only_limiting_partner"],
            f"A {meta['seconds_null_a']}s B {meta['seconds_null_b']}s")
        save(out, args)

    for scen in args.scenarios:
        rng = np.random.default_rng(SEEDS[scen])
        n = args.label_reps if scen == "labels" else args.reps
        reps = []
        for rep in range(n):
            if scen == "labels":
                ad, shuffled = shuffle_labels(adata, rng), None
            else:
                shuffled = rec_only if scen == "receptor" else lig_only
                ad = shuffle_genes(adata, shuffled, rng)
            with tempfile.TemporaryDirectory() as tmp:
                res, q, meta, _ = run(ad, "label", args.db, args.B, seed=rep + 1, workdir=Path(tmp))
            if scen == "labels":
                null = np.ones(len(res), bool)
            else:
                null = known_null(res, per, shuffled, f"{'receptor' if scen == 'receptor' else 'ligand'}_complex")
            tested = res.tested.to_numpy(bool)
            r = dict(meta, n_known_null=int(null.sum()), n_known_null_tested=int((null & tested).sum()), procs={})
            for m, v in q.items():
                r["procs"][m] = {str(al): dict(R=int((v <= al).sum()), V=int(((v <= al) & null).sum())) for al in ALPHAS}
            if rep == 0 and scen != "labels" and (null & tested).any():   # type-I error, p calibration
                rows = null & tested
                r["type1"] = {}
                for m, p in (("null_a", res.p_rand.fillna(1.0).to_numpy()), ("null_b", res.p_null_b.to_numpy())):
                    pk = p[rows]
                    r["type1"][m] = dict(rate={str(al): float((pk <= al).mean()) for al in ALPHAS},
                                         quantiles=np.quantile(pk, QS).tolist())
                lim = res.limiting_partner_null_b.to_numpy()
                fa = (q["null_a"] <= ALPHA) & null
                r["false_null_a_failing_partner_null_b"] = {k_: int((fa & (lim == k_)).sum()) for k_ in ("ligand", "receptor")}
            reps.append(r)
            log(args.tissue, scen, rep, {m: (d[str(ALPHA)]["R"], d[str(ALPHA)]["V"]) for m, d in r["procs"].items()},
                f"A {meta['seconds_null_a']}s B {meta['seconds_null_b']}s")
        summ = {}
        for m in reps[0]["procs"]:
            summ[m] = {}
            for al in map(str, ALPHAS):
                R = np.array([x["procs"][m][al]["R"] for x in reps])
                V = np.array([x["procs"][m][al]["V"] for x in reps])
                fdp = np.where(R > 0, V / np.maximum(R, 1), 0.0)
                if scen == "labels":
                    k = int((R > 0).sum())
                    summ[m][al] = dict(p_any_call=k / len(R), ci95=clopper_pearson(k, len(R)), mean_calls=float(R.mean()))
                else:
                    summ[m][al] = dict(fdr=mean_ci(fdp), mean_calls=float(R.mean()), mean_false=float(V.mean()))
        out["scenarios"][scen] = dict(summary=summ, reps=reps)
        log(args.tissue, scen, "summary at 0.05:", {m: v[str(ALPHA)] for m, v in summ.items()})
        save(out, args)
    return out


def part_sensitivity(args):
    adata = load_tissue(args.tissue, args.db)
    _, rec_only, per = side_genes(args.db)
    ad = shuffle_genes(adata, rec_only, np.random.default_rng(SEEDS["receptor"]))   # replicate 0 of calibration
    out = dict(tissue=args.tissue, scenario="receptor, replicate 0", by_B={})
    for B in args.Bs:
        with tempfile.TemporaryDirectory() as tmp:
            res, q, meta, _ = run(ad, "label", args.db, B, seed=1, workdir=Path(tmp))
        null = known_null(res, per, rec_only, "receptor_complex")
        out["by_B"][str(B)] = dict(meta, procs={m: {str(al): dict(R=int((v <= al).sum()), V=int(((v <= al) & null).sum()))
                                                    for al in ALPHAS} for m, v in q.items()})
        log("B", B, {m: d[str(ALPHA)] for m, d in out["by_B"][str(B)]["procs"].items()})
        save(out, args)
    return out


def fake_data(data: cpdb_null.CpdbData, X_all, genes_all, db_genes, rng):
    """CpdbData with the real interactions plus one fake copy of each simple-simple one: each partner's
    gene replaced by a non-CellPhoneDB gene from its mean-expression x detection stratum
    (ccc_nulls.nulls.gene_strata). Fake pairs are in no database: every fake call is false."""
    strata = gene_strata(X_all)
    pos = {g: i for i, g in enumerate(genes_all)}
    expressed = np.asarray((X_all > 0).sum(0)).ravel() > 0
    pool = expressed & ~pd.Index(genes_all).isin(list(db_genes))
    nrow = len(data.row_ids)
    inter = data.interactions
    simple = (data.i1 < nrow) & (data.i2 < nrow)
    new_rows, new_ids, fake = [], [], []
    for k in np.flatnonzero(simple):
        pair = []
        for f in (data.i1[k], data.i2[k]):
            g = data.row_genes[f]
            cand = np.flatnonzero(pool & (strata == strata[pos[g]])) if g in pos else np.flatnonzero(pool)
            new_rows.append(int(rng.choice(cand)))
            new_ids.append(f"fake:{len(new_ids)}")
            pair.append(new_ids[-1])
        fake.append(dict(multidata_1_id=pair[0], multidata_2_id=pair[1], id_cp_interaction=f"FAKE:{inter.id_cp_interaction.iloc[k]}"))
    Xf = sparse.csr_matrix(sparse.csc_matrix(X_all)[:, new_rows].T, dtype=np.float64)
    fd = cpdb_null.CpdbData(
        counts=sparse.vstack([data.counts, Xf], format="csr"), row_ids=np.concatenate([data.row_ids.astype(object), new_ids]),
        row_genes=np.concatenate([data.row_genes, np.asarray(genes_all)[new_rows]]), codes=data.codes,
        cell_types=data.cell_types, complexes=data.complexes,
        interactions=pd.concat([inter[["multidata_1_id", "multidata_2_id", "id_cp_interaction"]], pd.DataFrame(fake)],
                               ignore_index=True))
    return fd, int(simple.sum())


def part_tnbc(args):
    from ccc_nulls.evaluation import openproblems as opm
    from tools.cellchat.match_null_a import bootstrap_metrics

    name = "tnbc_wu2021"
    meta_ = datasets.DATASETS[name]
    groupby, keys, sample_key = meta_["groupby"], meta_["truth_keys"], meta_["sample_key"]
    full, truth = datasets.load_openproblems(name)
    full.obs[groupby] = full.obs[groupby].astype(str).str.replace("|", "/", regex=False).astype("category")
    db_genes = adapter.database_genes(args.db)
    adata = full[:, full.var_names.isin(db_genes)].copy()
    log(f"{name}: {adata.n_obs:,} cells, {adata.obs[groupby].nunique()} types, {adata.n_vars} CellPhoneDB genes")
    with tempfile.TemporaryDirectory() as tmp:
        res, q, meta, data = run(adata, groupby, args.db, args.B, seed=0, workdir=Path(tmp))
    out = dict(dataset=name, B=args.B, alpha=ALPHA, **meta, **describe(adata, res, q))
    log("calls:", out["calls"], out["overlap"], out["null_a_only_limiting_partner"])

    # Open Problems cytokine truth: each procedure gates CellPhoneDB's mean (as Open Problems does)
    rows = res[["ligand_complex", "target"]]
    score = res.score_obs.to_numpy()
    joined = {}
    for m in ["mean_ungated", "null_a", "null_b", "null_b|testable"]:
        pred = opm.predictions(rows, score, keys, p=None if m == "mean_ungated" else q[m], alpha=ALPHA)
        gt = opm._join(pred, truth, keys)
        joined[m] = (gt.response.to_numpy().astype(int), gt.score.to_numpy())
    rnd = truth[list(keys)].copy()
    rnd["score"] = np.random.default_rng(1).random(len(rnd))
    gt = opm._join(rnd, truth, keys)
    joined["random"] = (gt.response.to_numpy().astype(int), gt.score.to_numpy())
    out["openproblems"] = bootstrap_metrics(joined, ref="null_b", n_boot=args.metric_boot)
    for m, d in out["openproblems"]["methods"].items():
        log(f"  {m}: AUPRC {d['auprc']:.3f}, OR {d['odds_ratio']:.2f}",
            f"B - this {d['auprc_diff_ref_minus_this']:+.3f} {np.round(d['auprc_diff_ci95'], 3).tolist()}"
            if "auprc_diff_ci95" in d else "")
    save(out, args)

    # fake gene pairs (E2): one BH family over the tested real + fake rows
    fd, n_simple = fake_data(data, sparse.csr_matrix(full.X), full.var_names.astype(str).to_numpy(), db_genes,
                             np.random.default_rng(11))
    na, nb, _ = nulls_on(fd, args.B, seed=0)
    is_fake = np.char.startswith(na.interaction_ids.astype(str), "FAKE:")[:, None].repeat(len(na.pairs), 1).ravel()
    tested = na.tested.ravel()
    pa, pb = np.nan_to_num(na.p.ravel(), nan=1.0), nb["p"].ravel()
    out["fake_pairs"] = dict(n_fake_interactions=n_simple, n_tested_real=int((tested & ~is_fake).sum()),
                             n_tested_fake=int((tested & is_fake).sum()))
    for m, p in (("null_a", pa), ("null_b", pb)):
        c = over(p, tested) <= ALPHA
        out["fake_pairs"][m] = dict(real_calls=int((c & ~is_fake).sum()), fake_calls=int((c & is_fake).sum()),
                                    fake_rate_tested=float(c[tested & is_fake].mean()),
                                    real_rate_tested=float(c[tested & ~is_fake].mean()),
                                    fake_share_of_calls=float((c & is_fake).sum() / max(c.sum(), 1)))
    log("fake pairs:", out["fake_pairs"])
    del full
    save(out, args)

    # reproducibility across patient halves
    patients = np.array(sorted(adata.obs[sample_key].astype(str).unique()))
    half = set(np.random.default_rng(0).permutation(patients)[: len(patients) // 2])
    in_half = adata.obs[sample_key].astype(str).isin(half).to_numpy()
    sets = []
    for h, mask in enumerate([in_half, ~in_half]):
        sub = adata[mask].copy()
        vc = sub.obs[groupby].value_counts()
        sub = sub[sub.obs[groupby].isin(vc[vc >= 5].index)].copy()
        sub.obs[groupby] = sub.obs[groupby].astype(str).astype("category")
        log(f"half {h + 1}: {sub.n_obs:,} cells")
        with tempfile.TemporaryDirectory() as tmp:
            r, qh, _, _ = run(sub, groupby, args.db, args.B, seed=h + 1, workdir=Path(tmp))
        key = list(zip(r.native_interaction_id.astype(str), r.source.astype(str), r.target.astype(str)))
        sets.append({m: {k for k, c in zip(key, v <= ALPHA) if c} for m, v in qh.items()})
    out["halves"] = {m: dict(n1=len(sets[0][m]), n2=len(sets[1][m]),
                             jaccard=len(sets[0][m] & sets[1][m]) / max(len(sets[0][m] | sets[1][m]), 1)) for m in sets[0]}
    log("halves:", {m: round(v["jaccard"], 3) for m, v in out["halves"].items()})
    return out


def save(out, args):
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / f"{args.part}_{args.tissue if args.part != 'tnbc' else 'wu2021'}.json"
    f.write_text(json.dumps(out, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x)))
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="CellPhoneDB database zip (v5.0.0)")
    ap.add_argument("--part", required=True, choices=["calibration", "sensitivity", "tnbc"])
    ap.add_argument("--tissue", default="mammary")
    ap.add_argument("--scenarios", nargs="+", default=["receptor", "ligand", "labels"])
    ap.add_argument("--reps", type=int, default=20)
    ap.add_argument("--label-reps", type=int, default=10)
    ap.add_argument("--skip-real", action="store_true")
    ap.add_argument("--B", type=int, default=1000)
    ap.add_argument("--Bs", type=int, nargs="+", default=[500, 1000, 5000])
    ap.add_argument("--metric-boot", type=int, default=1000)
    args = ap.parse_args()
    t0 = time.time()
    out = {"calibration": part_calibration, "sensitivity": part_sensitivity, "tnbc": part_tnbc}[args.part](args)
    out["seconds"] = round(time.time() - t0, 1)
    log("wrote", save(out, args))


if __name__ == "__main__":
    main()
