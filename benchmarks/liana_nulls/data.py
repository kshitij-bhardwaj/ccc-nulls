"""The datasets of the README, loaded once and the same way for every tool.

    tnbc_wu2021             human, 42,512 cells, 29 types; truths: Visium colocalisation, cytokines
    allen_brain_atlas       mouse, 14,249 cells, 23 types; truth: Open Problems spatial source-target
    tabula_sapiens_mammary  human, calibration and power (no truth)
    tabula_sapiens_liver    human, calibration (no truth)
    moffitt2018_hypothalamus  mouse, scRNA of Moffitt et al. 2018; truth: MERFISH neighbourhoods (Amendment 1)

Cell-type names: "|" becomes "/" (CellPhoneDB's pair separator), in the data and in every truth.
"""
from __future__ import annotations

from ccc_nulls import datasets

SPECIES = {"tnbc_wu2021": "human", "allen_brain_atlas": "mouse",
           "tabula_sapiens_mammary": "human", "tabula_sapiens_liver": "human", "moffitt2018_hypothalamus": "mouse"}


def clean(labels):
    return labels.astype(str).str.replace("|", "/", regex=False).astype("category")


def database_genes(species: str = "human") -> set:
    """Genes any of the tools or LIANA+'s consensus can use (keeps Tabula Sapiens within memory)."""
    import liana as li

    from tools.cellchat import cellchat_py as cc
    from tools.cellphonedb import adapter as cpdb

    genes = set(cc.extract_gene(cc.subset_db(cc.load_db(species))))
    if species == "human":
        genes |= set(cpdb.database_genes(datasets.data_dir() / "cellphonedb_v5.0.0" / "cellphonedb.zip"))
    res = li.rs.select_resource("consensus" if species == "human" else "mouseconsensus")
    for col in ("ligand", "receptor"):
        genes |= {g for x in res[col] for g in str(x).split("_")}
    return genes


def load(name: str, extra_genes=()):
    """(adata with X = ln(1 + CP10k), groupby, extras) where extras holds the dataset's truths.
    Tabula Sapiens is cut to the database genes plus ``extra_genes`` (memory)."""
    extras = {"species": SPECIES[name]}
    if name in ("tnbc_wu2021", "allen_brain_atlas"):
        adata, truth = datasets.load_openproblems(name)
        groupby = datasets.DATASETS[name]["groupby"]
        adata.obs[groupby] = clean(adata.obs[groupby])
        for c in truth.columns.intersection(["source", "target"]):
            truth[c] = truth[c].astype(str).str.replace("|", "/", regex=False)
        extras["openproblems_truth"] = truth
        if name == "allen_brain_atlas":
            adata.uns["ccc_target"] = truth
    elif name == "moffitt2018_hypothalamus":
        from .hypothalamus import build
        adata, groupby = build(), "cell_class"
    else:
        adata = datasets.load_tabula_sapiens(name, genes=database_genes("human") | set(extra_genes))
        groupby = "label"
        adata.obs[groupby] = clean(adata.obs[groupby])
    return adata, groupby, extras


def orthologs(species: str = "mouse"):
    """One-to-one human <-> ``species`` gene symbols from HCOP (>= 3 supporting resources), cached in data/."""
    import pandas as pd

    f = datasets.data_dir() / f"hcop_human_{species}_1to1.csv"
    if not f.exists():
        import liana as li

        m = li.rs.get_hcop_orthologs(target_organism=species, min_evidence=3)
        m = m[["human_symbol", f"{species}_symbol"]].dropna().drop_duplicates()
        one = (~m.human_symbol.duplicated(keep=False)) & (~m[f"{species}_symbol"].duplicated(keep=False))
        m[one].rename(columns={f"{species}_symbol": "other_symbol"}).to_csv(f, index=False)
    return pd.read_csv(f)


def to_human(adata, species: str = "mouse"):
    """A copy with genes renamed to their 1:1 human orthologs; genes without one are dropped."""
    m = orthologs(species)
    mp = dict(zip(m.other_symbol, m.human_symbol))
    keep = [g for g in adata.var_names.astype(str) if g in mp]
    out = adata[:, keep].copy()
    out.var_names = [mp[g] for g in keep]
    return out, dict(source="HCOP, min_evidence = 3, one-to-one", genes_in=int(adata.n_vars), genes_mapped=len(keep))


def decoy_universe(name: str, n: int = 5000, seed: int = 0) -> set:
    """Amendment 3: for Tabula Sapiens (loaded as database genes only), ``n`` genes outside every
    database, drawn with ``seed``, to load as well so B2 has decoys. Empty for the other data sets,
    which are loaded with all their genes."""
    if not name.startswith("tabula_sapiens_"):
        return set()
    import h5py
    import numpy as np

    with h5py.File(datasets.path(name), "r") as f:
        symbols = np.unique(datasets._h5_strings(f, "raw/var/feature_name").astype(str))
    db = {g.upper() for g in database_genes("human") | database_genes("mouse")}
    pool = np.array([g for g in symbols if g.upper() not in db])
    rng = np.random.default_rng(seed)
    return set(rng.choice(pool, min(n, len(pool)), replace=False).tolist())
