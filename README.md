# ccc-nulls

**Structure-preserving null models and false-discovery control for ligand–receptor cell–cell communication (CCC) inference.**

Ligand–receptor tools such as CellPhoneDB, CellChat and LIANA usually assess
significance by shuffling cell-type labels uniformly. That null tests whether cell type
affects a gene's expression *at all*, and for nearly every gene it does. As a result, abundant
ligands get declared "significant" whether or not anything specific is going on.

`ccc-nulls` is a **tool-agnostic significance layer**. It takes any tool's ranked LR table
plus the expression matrix it was computed from, and returns calibrated p- and q-values under two nulls:

| | Null A: random | Null B: structure-preserving |
|---|---|---|
| What is permuted | Cell-type labels, uniformly | Gene identity, within mean-expression × detection-rate strata |
| What is preserved | Cell-type sizes | Cell-type sizes, every gene's profile across cell types, and abundance |
| Question answered | "Does cell type matter for these genes?" | "Is this pair extreme *for genes of its abundance*?" |
| Used by | CellPhoneDB, CellChat (default) | This package |

It composes with existing tools instead of replacing them. The LR scoring stays the tool's own.

---

## Installation

```bash
pip install "git+https://github.com/kshitij-bhardwaj/ccc-nulls.git"
# with the liana-py / scanpy integration:
pip install "ccc-nulls[liana] @ git+https://github.com/kshitij-bhardwaj/ccc-nulls.git"
```

For development:

```bash
git clone https://github.com/kshitij-bhardwaj/ccc-nulls.git && cd ccc-nulls
pip install -e ".[test]" && pytest
```

Requires Python ≥ 3.10, with numpy, scipy and pandas. anndata, scanpy and liana are optional.

## Quickstart

```python
import ccc_nulls as cn

# candidates: DataFrame with columns ligand, receptor, source, target (anything else is kept)
res = cn.significance_layer(
    X, genes, labels, candidates,      # X: cells x genes, normalised (not scaled)
    score="mean",                      # "mean" (CellPhoneDB) | "product" | callable
    B=1000,                            # Null A permutations; Null B is exact
    resource_genes=lr_resource_genes,  # kept out of Null B's decoy pool
)
cn.overlap_summary(res)                # both / A_only / B_only / none at q < 0.05
```

### With liana-py

```python
import liana as li, ccc_nulls as cn

li.mt.cellphonedb(adata, groupby="cell_type", use_raw=True, key_added="cpdb")
cands = cn.adapters.from_liana(adata.uns["cpdb"], method="cellphonedb")
res = cn.significance_layer_anndata(adata, "cell_type", cands, use_raw=True)
```

Columns added to your table:

| Column | Meaning |
|---|---|
| `score_obs` | The tool score recomputed from `X` (both nulls test this statistic) |
| `p_rand`, `q_rand` | Null A p-value and BH q-value |
| `z_strat`, `p_strat`, `q_strat` | Null-A-standardised score, Null B p-value, BH q-value |
| `p_strat_min` | Smallest p Null B can reach for this row (its resolution) |
| `flag` | `both`, `A_only`, `B_only`, `none`. **`A_only` = significant only because the uniform null ignored abundance** |

See [`examples/quickstart.py`](examples/quickstart.py) (simulation with known truth) and
[`examples/liana_pbmc.py`](examples/liana_pbmc.py) (liana-py end to end, no download).

## Method

For a candidate (L, R, X, Y) with tool score `s`:

**Null A.** Draw B uniform permutations π_b of the cell labels and compute
`p_rand = (1 + #{b : s(π_b) ≥ s_obs}) / (B + 1)`.

**Null B.**
1. **Strata.** Bin genes by mean-expression decile, then by detection-rate tercile within each decile (30 strata).
2. **Standardise.** Centre and scale the tool score by its *exact* Null-A moments. A permuted
   cluster is a sample drawn without replacement, so `Var = σ²/n_k · (n − n_k)/(n − 1)`. This
   puts genes of every abundance on one scale at no permutation cost.
3. **Decoys.** For L, take background genes (not in the LR resource) from L's stratum that pass
   the same expression filter in X as the candidate did. Do the same for R in Y.
4. **Exact p-value.** Under H0, (L, R) is a uniform draw from (decoys_L ∪ {L}) × (decoys_R ∪ {R}).
   `p_strat` is the candidate's upper-tail rank in that product set. The pairs that reuse L or R
   are included, which keeps the test valid when one gene is extreme.

**Corroboration (planned).** `cauchy_combine` (Liu & Xie, 2020) merges `p_strat` with a
NicheNet-style target-enrichment p-value. It is valid under the dependence between the two.

## Benchmarks

`benchmarks/simulation_study.py` uses negative-binomial data with 8,000 genes, 6 cell types,
400 LR pairs and about 9,400 tests per run. Results below are means over 5 seeds, with Null A
at B = 10,000, from [`benchmarks/results/simulation.json`](benchmarks/results/simulation.json).

| Regime | Ground truth | Null A | Null B |
|---|---|---|---|
| S0: fully permuted | No structure, no CCC | FPR 5.3% | FPR 4.8% |
| S1: structured null | Cell identity shifts every gene (τ = 0.25); no CCC | **FPR 29.7%**, ~2,200 false calls/run | **FPR 5.2%**, 0 calls |
| S2: planted, LFC 1.0 | 60 true events | FDR 97%, power 100% | FDR 6%, power 74% |
| S2: planted, LFC 1.5 | | FDR 97%, power 100% | FDR 9%, power 90% |
| S2: planted, LFC 0.5 | | FDR 98%, power 93% | FDR 0%, power 0% |

Same picture with the product score (S1: 24.8% vs 5.1%). Null A's false calls concentrate in
the most abundant ligands: 3.2× more in the top three deciles than in the bottom three.

`benchmarks/binned_label_check.py` checks two common shortcuts:
- **Permuting cell labels within expression bins is not a fix.** Each gene still sees a uniform
  shuffle, so its null equals Null A's (log-p correlation 0.997, same 29% FPR).
- **Fully permuted data (S0) cannot validate a null.** Both nulls look calibrated there. Use a
  structured null.

## Limitations

- **Power:** Null B gives up weak effects (LFC 0.5) in exchange for error control. Report both q-values.
- **Resolution:** Null B's smallest p is `1 / ((|decoys_L|+1)(|decoys_R|+1))`. Run it on the full
  transcriptome. On highly-variable-gene subsets the strata are too small (watch `p_strat_min`).
  The bundled PBMC example shows this.
- **Assumption:** the S1 structure level τ = 0.25 is assumed. Estimating it on real atlases is ongoing.
- **Scores:** standardisation supports the built-in `mean` and `product` scores. For custom scores
  pass `standardize=False` through the low-level `null_b`.
- **Adapters:** liana-py output is supported now (it covers CellPhoneDB, CellChat, NATMI,
  Connectome, SingleCellSignalR, logFC and rank aggregate). Native CellPhoneDB/CellChat output is on the roadmap.

## Repository layout

```
src/ccc_nulls/
  core.py          significance_layer, overlap_summary, flag_overlap, aggregate_complexes
  nulls.py         null_a, null_b, gene_strata, perm_moments, null_b_binned_labels (diagnostic)
  scores.py        CellPhoneDB mean, product, standardisation
  stats.py         Benjamini–Hochberg, Cauchy combination
  adapters.py      from_liana, to_schema
  io.py            AnnData helpers
  simulate.py      simulator with known truth (S0 / S1 / S2)
benchmarks/        simulation_study.py, binned_label_check.py, results/*.json
examples/          quickstart.py, liana_pbmc.py
tests/             pytest suite (calibration, exactness, API)
```

## References

1. Dimitrov D. et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat Commun* 13, 3224 (2022).
2. Dimitrov D. et al. LIANA+ provides an all-in-one framework for cell–cell communication inference. *Nat Cell Biol* 26, 1613–1622 (2024).
3. Open Problems in Single-Cell Analysis: cell–cell communication benchmark.
4. Jin S. et al. Inference and analysis of cell-cell communication using CellChat. *Nat Commun* 12, 1088 (2021).
5. Efremova M. et al. CellPhoneDB. *Nat Protoc* 15, 1484–1506 (2020); Troulé K. et al. CellPhoneDB v5. *Nat Protoc* (2025).
6. Browaeys R., Saelens W., Saeys Y. NicheNet. *Nat Methods* 17, 159–162 (2020).
7. Liu Y., Xie J. Cauchy combination test. *JASA* 115, 393–402 (2020).

## Authors

Kaushal, Kshitij Bhardwaj, Ankit, Aritra. Developed as CS690A (Computational Genomics), Project 6.
