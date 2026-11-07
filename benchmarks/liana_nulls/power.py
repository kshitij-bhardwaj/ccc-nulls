"""T2: planted edges on a background where no edge is true. Power and FDR of Null A and Null B.

Background: Tabula Sapiens Mammary, 10,000 random cells, cell-type labels shuffled. Each replicate
plants 30 L-R pairs, single-gene on both sides and in both tools' databases (CellChatDB and
CellPhoneDB), with disjoint genes. Each gets a random (s, t), s != t: the ligand is raised in s and
the receptor in t, by fold change f on the CP10k scale (non-zero x -> log1p(f expm1(x)); a zero is
switched on with probability 1 - 1/f, drawing a non-zero value of the gene in that type, raised the
same way). Every effect size uses the same cells, labels and pairs (paired across f).

Truth of a row (s, t, L, R) for the claim "s signals to t":
- true: every ligand subunit raised in s and every receptor subunit raised in t;
- one-sided null: one side raised in its type, the other not (the false call Null A makes);
- ambiguous: a complex partner with some but not all subunits raised; left out of FDP and power;
- null: everything else.
Power = true calls / true (s, t, L, R) tuples planted in the data (a tool that filters a true row out
misses it). FDP = calls on non-true rows / calls, ambiguous rows excluded from both.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

N_CELLS = 10_000
N_PAIRS = 30
EFFECTS = (1.5, 2.0, 3.0, 4.0)


def candidate_pairs(adata, cpdb_db, min_frac: float = 0.02) -> pd.DataFrame:
    """Single-gene L-R pairs in both CellChatDB and CellPhoneDB, both genes detected in >= min_frac of cells."""
    from tools.cellchat import cellchat_py as cc
    from tools.cellphonedb import adapter as cpdb

    inter = cc.subset_db(cc.load_db("human"))["interaction"]
    genes = set(adata.var_names.astype(str))
    cc_pairs = inter[inter.ligand.isin(genes) & inter.receptor.isin(genes)
                     & (inter.agonist == "") & (inter.antagonist == "")
                     & (inter.co_A_receptor == "") & (inter.co_I_receptor == "")][["ligand", "receptor"]]
    pairs = cc_pairs.merge(cpdb.resource_pairs(cpdb_db), on=["ligand", "receptor"]).drop_duplicates()
    X = sparse.csc_matrix(adata.X)
    frac = pd.Series(np.diff(X.indptr) / X.shape[0], index=adata.var_names.astype(str))
    ok = (frac.reindex(pairs.ligand).to_numpy() >= min_frac) & (frac.reindex(pairs.receptor).to_numpy() >= min_frac)
    return pairs[ok].reset_index(drop=True)


def pick_pairs(cands: pd.DataFrame, n: int, rng) -> pd.DataFrame:
    """n pairs with disjoint genes, in random order."""
    used, out = set(), []
    for i in rng.permutation(len(cands)):
        l, r = cands.ligand.iat[i], cands.receptor.iat[i]
        if l in used or r in used or l == r:
            continue
        used |= {l, r}
        out.append((l, r))
        if len(out) == n:
            break
    return pd.DataFrame(out, columns=["ligand", "receptor"])


def background(adata, groupby: str, rep: int):
    """10,000 random cells with their labels shuffled: no cell-type structure is left."""
    rng = np.random.default_rng(1000 + rep)
    keep = np.sort(rng.choice(adata.n_obs, min(N_CELLS, adata.n_obs), replace=False))
    bg = adata[keep].copy()
    bg.obs[groupby] = pd.Categorical(rng.permutation(bg.obs[groupby].astype(str).to_numpy()))
    vc = bg.obs[groupby].value_counts()
    return bg[bg.obs[groupby].isin(vc[vc >= 20].index)].copy()


def design(bg, groupby: str, cands: pd.DataFrame, rep: int) -> pd.DataFrame:
    rng = np.random.default_rng(2000 + rep)
    pairs = pick_pairs(cands, N_PAIRS, rng)
    types = np.array(sorted(bg.obs[groupby].astype(str).unique()))
    st = np.array([rng.choice(types, 2, replace=False) for _ in range(len(pairs))])
    return pairs.assign(source=st[:, 0], target=st[:, 1])


def raise_gene(X: sparse.csc_matrix, j: int, cells: np.ndarray, f: float, rng) -> None:
    """In place: column j raised by fold change f in ``cells`` (log1p CP10k)."""
    col = np.asarray(X[:, j].todense()).ravel()
    x = col[cells]
    nz = x > 0
    pool = x[nz] if nz.any() else col[col > 0]
    on = ~nz & (rng.random(len(x)) < 1 - 1 / f)
    x[on] = rng.choice(pool, on.sum()) if len(pool) else 0.0
    x[x > 0] = np.log1p(f * np.expm1(x[x > 0]))
    col[cells] = x
    X[:, j] = sparse.csc_matrix(col[:, None])


def plant(bg, groupby: str, plan: pd.DataFrame, f: float, rep: int):
    """A copy of ``bg`` with each planned pair's ligand raised in its source and receptor in its target."""
    rng = np.random.default_rng(3000 + rep)
    a = bg.copy()
    X = sparse.csc_matrix(a.X, dtype=np.float32)
    var = pd.Index(a.var_names.astype(str))
    labels = a.obs[groupby].astype(str).to_numpy()
    for p in plan.itertuples():
        raise_gene(X, var.get_loc(p.ligand), np.where(labels == p.source)[0], f, rng)
        raise_gene(X, var.get_loc(p.receptor), np.where(labels == p.target)[0], f, rng)
    a.X = sparse.csr_matrix(X)
    return a


def label_truth(rows: pd.DataFrame, plan: pd.DataFrame) -> np.ndarray:
    """'true' / 'one_sided' / 'ambiguous' / 'null' for every row (see module docstring)."""
    raised = {}                      # gene -> cell type it was raised in
    for p in plan.itertuples():
        raised[p.ligand], raised[p.receptor] = p.source, p.target

    def side(complex_names, types):
        out = np.empty(len(complex_names), dtype=object)
        for i, (cx, ty) in enumerate(zip(complex_names, types)):
            sub = str(cx).split("_")
            hit = [raised.get(g) == ty for g in sub]
            out[i] = "yes" if all(hit) else ("partial" if any(hit) else "no")
        return out

    L = side(rows.ligand_complex.to_numpy(), rows.source.astype(str).to_numpy())
    R = side(rows.receptor_complex.to_numpy(), rows.target.astype(str).to_numpy())
    lab = np.full(len(rows), "null", dtype=object)
    lab[(L == "yes") ^ (R == "yes")] = "one_sided"
    lab[(L == "partial") | (R == "partial")] = "ambiguous"
    lab[(L == "yes") & (R == "yes")] = "true"
    return lab


def metrics(rows: pd.DataFrame, truth: np.ndarray, is_planted: np.ndarray, n_planted: int, call_cols: dict) -> dict:
    """Per call rule: power on the planted tuples, FDP over non-ambiguous rows (true = both sides raised)."""
    keep = truth != "ambiguous"
    out = {}
    for name, col in call_cols.items():
        c = rows[col].to_numpy(bool) & keep
        R = int(c.sum())
        false = c & (truth != "true")
        out[name] = dict(calls=R, planted_calls=int((c & is_planted).sum()), power=float((c & is_planted).sum() / n_planted),
                         fdp=float(false.sum() / R) if R else 0.0,
                         one_sided_false_calls=int((c & (truth == "one_sided")).sum()),
                         null_false_calls=int((c & (truth == "null")).sum()))
    out.update(n_rows=int(keep.sum()), n_true_rows=int((truth == "true").sum()), n_planted_rows=int(is_planted.sum()),
               n_one_sided_rows=int((truth == "one_sided").sum()), n_ambiguous_rows=int((~keep).sum()))
    return out
