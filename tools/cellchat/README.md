# CellChat

One owner, assigned within the team. Only the owner edits this folder; anything shared goes in
`src/ccc_nulls/` through a reviewed pull request.

| | |
|---|---|
| Tool | CellChat **2.1.2**, the latest release (`jinworks/CellChat`, tag `v2.1.2`, April 2025), **ported to Python**: no R needed |
| Also available | `2.2.0.9001` (GitHub main, `75253cd`) and `1.0.0` (the 2021 paper, `sqjin/CellChat`, tag `v1.0.0`) |
| Statistic tested by the nulls | `prob`, CellChat's communication probability (LIANA's re-implementation calls it `lr_probs`) |
| Native null (Null A) | Yes: `computeCommunProb` permutes the cell labels, `nboot = 100`, `seed.use = 1` |
| Extra Python dependency | `rdata` (reads CellChatDB from the package's `.rda` once): `pip install rdata` |

## Files

| File | What it is |
|---|---|
| `cellchat_py.py` | CellChat ported line by line from R: CellChatDB, `subsetDB`, `subsetData`, `identifyOverExpressedGenes` (presto's Wilcoxon test), `identifyOverExpressedInteractions`, `computeCommunProb` (2.1.2, 2.2.0.9001, 1.0.0), and R's random-number generator (`set.seed`, `sample.int`) |
| `adapter.py` | dataset → CellChat inputs (`prepare_inputs`) → CellChat run (`run_native`, the port) → common schema (`native_to_schema`) |
| `nulls.py` | `null_a`: CellChat's permutation test, vectorised; `null_b`: to do |
| `match_null_a.py` | TNBC: port vs Null A for every version, independent permutations, version differences, ground-truth metrics with bootstrap CIs |
| `calibrate_null_a.py` | Tabula Sapiens (Mammary, Liver): full-permutation calibration |
| `results/` | the JSON outputs of the two scripts |

## Versions

All three use the same CellChatDB file (identical in v2.1.2 and 75253cd) and the same filters.

| | 1.0.0 (paper) | **2.1.2 (release, default)** | 2.2.0.9001 (main) |
|---|---|---|---|
| complexes, co-receptors | per cell, then trimean per group | trimean per gene, then combined | as 2.1.2 |
| quantile in `triMean` | `stats::quantile`: `(1 - h) a + h b`, only where `b != a` | `stats::quantile` | `collapse::fquantile`: `a + h (b - a)` |
| p-value | `#(prob_perm >= prob) / nboot` | `#(prob_perm > prob) / nboot` | as 2.1.2 |
| default data | smoothed over a protein network (`projectData`); **we run it on raw data** | raw | raw |

2.1.2 and 2.2.0.9001 differ only in the last bit of some quantiles: on TNBC their probabilities
differ in 1,186 of 1,757,690 rows by at most 6 × 10⁻¹⁷, and every p-value is the same.

## What CellChat computes, and what Null A re-implements

From `computeCommunProb` (R/modeling.R, v2.1.2), on `data = data.signaling / max(data.signaling)`:

| Step | CellChat | ours |
|---|---|---|
| group average | `triMean`: `mean(stats::quantile(x, c(.25, .5, .5, .75)))` per gene and group | same order statistics and interpolation, R's two-pass `mean` |
| ligand / receptor | single gene: its average; complex: geometric mean of the subunits' averages | same |
| co-receptors | `R * prod(1 + avg[co-A]) / prod(1 + avg[co-I])` | same, sequential products |
| law of mass action + Hill | `P1 = LR / (Kh + LR)`, `LR[s, t] = L[s] R[t]`, `Kh = 0.5`, `n = 1` | same |
| agonist / antagonist | `P2 = ag agᵀ`, `ag = prod(1 + avg / (Kh + avg))`; `P3 = an anᵀ`, `an = prod(Kh / (Kh + avg))` | same |
| population size | `P4`, off by default (`population.size = FALSE`) | `population_size=False` |
| **Null A** | `set.seed(seed.use)`, `nboot` permutations of all cells' labels (`sample.int`), recompute `prob`; `p = #(prob_perm > prob) / nboot`; `p = 1` where `prob = 0` | same; `rng="cellchat"` draws CellChat's own permutations with R's generator |
| significance | `subsetCommunication(thresh = 0.05)`: `p < 0.05`. **No q-values, no FDR control** | BH q-values added, for the comparison |

- CellChat's p has no `+1`, so it can be exactly 0, and with `nboot = 100` its resolution is
  0.01. `null_a` returns `p_null_a_cellchat` (CellChat's formula) and `p_null_a`
  (`(1 + #) / (B + 1)`, the repository's convention).
- Before Null A, CellChat keeps genes over-expressed in at least one group (Wilcoxon p < 0.05,
  logFC > 0) and the L-R pairs whose ligand and receptor both pass.
- LIANA's `cellchat` is a simplified re-implementation (no co-receptors, agonists or
  antagonists; LIANA's resource), so it is a reference, not the original.

## How the match is checked

```bash
python -m tools.cellchat.match_null_a --quick     # 8,000 TNBC cells
python -m tools.cellchat.match_null_a             # all 42,512 TNBC cells (~15 min)
pytest tools/cellchat                             # includes exact-match tests on 700 PBMCs
```

1. **R's generator.** `RRandom` reproduces known R outputs (`set.seed(1); sample(10)` is
   `9 4 7 1 2 5 3 10 6 8`, and four more), so the permutations are the ones CellChat draws.
2. **Reference.** `cellchat_py.compute_commun_prob` follows the R loop literally (dense sort per
   group, one L-R pair and one permutation at a time).
3. **Our Null A.** `nulls.null_a(..., B=100, seed=1)` is a separate, vectorised implementation
   (sparse order statistics, all pairs at once). It must match the reference exactly.
4. **Same null distribution.** `null_a(..., rng="numpy")` with independent permutations must
   agree with CellChat up to Monte-Carlo error.

### How far the match goes

| Claim | Established? | Evidence |
|---|---|---|
| our fast Null A = our line-by-line port | **yes, exactly** | prob, p, BH q and calls identical on every TNBC row, for all three versions |
| our permutations = the ones R draws for `set.seed(1)` | **yes** | `RRandom` reproduces 5 known R outputs |
| the port = CellChat's R source | **by reading the code** | every formula and operation order follows R/modeling.R of the version named |
| our numbers = what CellChat prints when run in R | **not checked** | CellChat in R was never run (the project is Python only) |

What could still differ from running CellChat in R, and why it should not:
- `log` / `exp`: numpy matches the C library bit for bit (checked on 10⁶ values each).
- R's `mean` uses long double; on Apple Silicon long double is double, which the port assumes.
  On x86 Linux R's `mean` may differ in the last bit, which can only flip an exact tie between a
  permuted and an observed probability.
- presto's p-values at exactly the 0.05 boundary can only change which L-R pairs enter, not Null A.

## Results: Null A on TNBC (Open Problems `tnbc_wu2021`, 42,512 cells, 29 groups)

`results/null_a_match_tnbc_wu2021_full.json` (`_quick` for 8,000 cells). CellChat's pipeline:
1,147 CellChatDB genes in the data, 1,143 over-expressed, 2,090 L-R pairs, 29 × 29 × 2,090 =
1,757,690 (pair, source, target) rows. The trimean is 0 unless more than ~25% of a group
expresses the gene, so most rows are 0 and get p = 1.

| Check | 2.1.2 (release) | 2.2.0.9001 | 1.0.0 |
|---|---|---|---|
| rows with prob > 0 | 14,872 | 14,872 | 14,031 |
| prob identical, port vs Null A (all rows) | **1,757,690 / 1,757,690** | 1,757,690 / 1,757,690 | 1,757,690 / 1,757,690 |
| p identical (rows with prob > 0) | **14,872 / 14,872** | 14,872 / 14,872 | 14,031 / 14,031 |
| BH q identical | **14,872 / 14,872** | 14,872 / 14,872 | 14,031 / 14,031 |
| calls at p < 0.05, ours / CellChat's | 13,024 / 13,024, 0 differ | 13,024, 0 differ | 12,470, 0 differ |
| p = 0 (no permutation reached the observed prob) | 12,862 | 12,862 | 12,300 |

Release 2.1.2 with independent permutations (same null, different draws):

| | within 2 Monte-Carlo sd | within 3 sd | calls at p < 0.05 (both / only ours / only CellChat) |
|---|---|---|---|
| B = 100 | 99.8% | 99.99% | 13,002 / 17 / 22 |
| B = 1,000 | 99.6% | 100% | 13,007 / 18 / 17 |

1.0.0 vs 2.1.2: 3,260 of the 1,757,690 probabilities differ, all involving complexes or
co-receptors; 13,765 of the 14,031 shared p-values are identical and 13,830 are ≥ 2.1.2's (ties
count against in 1.0.0); calls at p < 0.05: 12,371 in both, 653 only in 2.1.2, 99 only in 1.0.0.

CellChat's own significance on TNBC: 13,024 of 14,872 scored interactions (88%) have p < 0.05,
12,862 have p = 0, so their BH q-values are 0 as well and BH removes nothing. CellChat reports
no FDR, and its p-value cannot be read as one.

## Evaluation against ground truth (plan step 1, TNBC cytokine truth)

Field `openproblems` of the same file. Open Problems ligand-target protocol (812 truth rows, 166
positive); point estimates equal `ccc_nulls.evaluation.openproblems`; 95% CIs from 1,000
bootstrap resamples of the truth rows, paired across methods.

| | truth rows scored > 0 | AUPRC [95% CI] | odds ratio, top 5% [95% CI] | AUPRC, 2.1.2 minus this [95% CI] |
|---|---|---|---|---|
| **CellChat 2.1.2** | 41 | **0.350** [0.258, 0.447] | 3.45 [1.73, 6.43] | |
| CellChat 2.1.2, score kept where Null A p ≤ 0.05 | 41 | 0.350 [0.258, 0.447] | 3.45 | 0 |
| CellChat 2.2.0.9001 | 41 | 0.350 | 3.45 | 0 |
| CellChat 1.0.0 (paper) | 25 | 0.291 [0.193, 0.404] | 1.95 [0.77, 3.31] | 0.059 [0.001, 0.125] |
| LIANA's CellChat re-implementation | 158 | 0.235 [0.189, 0.293] | 3.09 [1.32, 5.54] | 0.115 [0.039, 0.191] |
| random | 812 | 0.214 [0.175, 0.264] | 1.51 [0.57, 3.08] | 0.136 [0.040, 0.231] |

- CellChat 2.1.2 ranks the cytokine truth better than LIANA's re-implementation and than random
  (the paired CIs exclude 0), and slightly better than 1.0.0 (CI just above 0).
- **Caveat:** CellChat gives a score > 0 to only 41 of the 812 truth rows (1.0.0: 25); the other
  ~770 are tied at the bottom. Its AUPRC rests on those few rows, hence the wide CI, and the top
  5% (40 rows) is essentially its whole scored set.
- Gating by Null A changes nothing: every top-scoring interaction already has p < 0.05.

## Calibration (plan step 2, Tabula Sapiens, full permutation)

`python -m tools.cellchat.calibrate_null_a --tissue all [--version 1.0.0]` →
`results/calibration_<tissue>_<version>.json`. Tabula Sapiens (the root README's calibration
collection), raw counts → ln(1 + CP10k), cell types with ≥ 5 cells:
**Mammary** 30,934 cells, 16 types, 2 donors; **Liver** 22,214 cells, 24 types, 4 donors.
Every cell's label is shuffled once (20 independent shuffles per tissue), so no interaction is
real, and CellChat runs end to end with its defaults (nboot = 100, seed.use = 1). Rates are
pooled over shuffles; 95% CIs from a bootstrap over shuffles.

P(p < 0.05) on the rows CellChat reports (prob > 0); a calibrated test gives 5/101 = **4.95%**:

| Fully shuffled labels | Mammary, 2.1.2 | Liver, 2.1.2 | Mammary, 1.0.0 | Liver, 1.0.0 |
|---|---|---|---|---|
| CellChat pipeline | **15.1%** [13.7, 16.4] | **12.6%** [10.9, 14.3] | **15.7%** [14.3, 17.2] | **12.9%** [11.0, 14.9] |
| without the over-expression filter | 9.4% [8.4, 10.3] | 9.4% [8.6, 10.3] | 9.8% [8.8, 10.9] | 9.7% [8.9, 10.6] |
| rows where no permutation reached prob = 0 (no filter, 3 shuffles) | **4.92%** | **5.10%** | | |
| P(p = 0), pipeline (expected 0.99%) | 4.3% [3.7, 4.8] | 3.2% [2.6, 3.8] | 4.5% [3.8, 5.2] | 3.2% [2.6, 4.0] |
| P(p < 0.05) over all rows (prob = 0 rows have p = 1) | 0.47% | 0.16% | 0.42% | 0.13% |
| BH false calls over 20 shuffles: CellChat's p / (1 + #)/(B + 1) | 3,667 / 0 | 4,157 / 0 | 3,283 / 0 | 3,524 / 0 |

With the real labels, 79% (Mammary) and 84% (Liver) of scored rows have p < 0.05.

What this shows:
1. **Our Null A is correct.** Where nothing is selected it is calibrated (4.92% and 5.10% vs
   4.95%), and over all rows it is conservative, as a permutation test must be.
2. **CellChat's reported p-values are anti-conservative, 2.5 to 3.2 times nominal, in both tissues
   and both versions**, because rows are selected with the same labels that are tested:
   - the over-expression filter keeps genes that look high in some group under these labels
     (15.1% → 9.4% on Mammary without it);
   - only rows with prob > 0 are reported, and a permutation that drives prob to 0 always loses
     to an observed prob > 0 (11.5% on Mammary where that happens, 4.9% where it does not).
3. **CellChat's p = 0 breaks BH.** With p = 0 allowed, BH calls 164 to 208 interactions per
   shuffled data set (mean over 20 shuffles) with no signal at all; with `(1 + #) / (B + 1)` it calls none.

Full permutation separates a valid test from an invalid one here only because CellChat's
selection is label-dependent; for abundance confounding, the fake-L-R-resource test is still
needed (docs/PLAN.md).

## Checklist

- [x] cite one exact CellChat version: release 2.1.2 (plus main and the paper's version)
- [x] keep all rows (prob = 0 rows get p = 1), not only significant ones, so nulls can be compared
- [x] complexes keep CellChat's names in `ligand` / `receptor`; genes are in `ligand_complex` / `receptor_complex`
- [x] `nulls.null_a`: reproduce CellChat's own permutation test; compare with `p_native`
- [x] calibration: 20 shuffles, two tissues, CIs; ground truth with bootstrap CIs
- [ ] `nulls.null_b`

## Run the tests and the checks

```bash
pytest tools/cellchat
python -m tools.cellchat.match_null_a --quick
python -m tools.cellchat.calibrate_null_a --download --tissue all
python -m tools.cellchat.calibrate_null_a --tissue all
```
