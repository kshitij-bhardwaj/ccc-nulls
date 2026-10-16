# CellChat

One owner, assigned within the team. Only the owner edits this folder; anything shared goes in
`src/ccc_nulls/` through a reviewed pull request.

| | |
|---|---|
| Tool | CellChat **2.2.0.9001**, `jinworks/CellChat@75253cd` (latest on GitHub, March 2026), **ported to Python**: no R needed |
| Statistic tested by the nulls | `prob`, CellChat's communication probability (LIANA's re-implementation calls it `lr_probs`) |
| Native null (Null A) | Yes: `computeCommunProb` permutes the cell labels, `nboot = 100`, `seed.use = 1` |
| Extra Python dependency | `rdata` (reads CellChatDB from the package's `.rda` once): `pip install -e ".[cellchat]"` |

## Files

| File | What it is |
|---|---|
| `cellchat_py.py` | CellChat's pipeline ported line by line from R: CellChatDB, `subsetDB`, `subsetData`, `identifyOverExpressedGenes` (presto's Wilcoxon test), `identifyOverExpressedInteractions`, `computeCommunProb`, and R's random-number generator (`set.seed`, `sample.int`) |
| `adapter.py` | dataset → CellChat inputs (`prepare_inputs`) → CellChat run (`run_native`, the port) → common schema (`native_to_schema`) |
| `nulls.py` | `null_a`: CellChat's permutation test, vectorised (fast); `null_b`: to do |
| `match_null_a.py` | runs the port and our Null A on a dataset and checks they agree |
| `results/` | match results (JSON) |

## What CellChat computes, and what Null A re-implements

From `computeCommunProb` (R/modeling.R at 75253cd), on `data = data.signaling / max(data.signaling)`:

| Step | CellChat | ours |
|---|---|---|
| group average | `triMean`: `mean(collapse::fquantile(x, c(.25, .5, .5, .75)))` per gene and group | same order statistics, collapse's `a + h (b - a)` interpolation, R's two-pass `mean` |
| ligand / receptor | single gene: its average; complex: geometric mean of the subunits' averages | same |
| co-receptors | `R * prod(1 + avg[co-A]) / prod(1 + avg[co-I])` | same, sequential products |
| law of mass action + Hill | `P1 = LR / (Kh + LR)`, `LR[s, t] = L[s] R[t]`, `Kh = 0.5`, `n = 1` | same |
| agonist / antagonist | `P2 = ag agᵀ`, `ag = prod(1 + avg / (Kh + avg))`; `P3 = an anᵀ`, `an = prod(Kh / (Kh + avg))` | same |
| population size | `P4`, off by default (`population.size = FALSE`) | `population_size=False` |
| **Null A** | `set.seed(seed.use)`, `nboot` permutations of all cells' labels (`sample.int`), recompute `prob`; `p = #(prob_perm > prob) / nboot`; `p = 1` where `prob = 0` | same; `rng="cellchat"` draws CellChat's own permutations with R's generator |
| significance | `subsetCommunication(thresh = 0.05)`: `p < 0.05`. **No q-values, no FDR control** | BH q-values added to both, for the comparison |

Notes:
- CellChat's p has no `+1`, so it can be exactly 0. With the default `nboot = 100` its resolution
  is 0.01. `null_a` returns both `p_null_a_cellchat` (CellChat's formula) and `p_null_a`
  (`(1 + #) / (B + 1)`, the repository's convention).
- CellChat filters candidates before Null A: genes over-expressed in at least one group
  (Wilcoxon p < 0.05, logFC > 0), and L-R pairs whose ligand and receptor both pass.
- LIANA's `cellchat` method is a simplified re-implementation (no co-receptors, agonists or
  antagonists, LIANA's own resource), so it is a reference, not the original.

## How the match is checked

```bash
python -m tools.cellchat.match_null_a --quick     # 8,000 TNBC cells (~3 min)
python -m tools.cellchat.match_null_a             # all 42,512 TNBC cells
pytest tools/cellchat                             # includes the exact-match test on 700 PBMCs
```

1. **R's generator.** `RRandom` reproduces known R outputs (`set.seed(1); sample(10)` is
   `9 4 7 1 2 5 3 10 6 8`, and four more), so the permutations are the ones CellChat draws.
2. **Reference.** `cellchat_py.compute_commun_prob` follows the R loop literally (dense sort per
   group, one L-R pair and one permutation at a time) and gives the prob and p-values CellChat
   reports.
3. **Our Null A.** `nulls.null_a(..., B=100, seed=1)` is a separate, vectorised implementation
   (sparse order statistics, all pairs at once). It must match the reference exactly: prob,
   p-values, BH q-values and calls.
4. **Same null distribution.** `nulls.null_a(..., rng="numpy")` with independent permutations
   (B = 100, 1,000) must agree with CellChat up to Monte-Carlo error.

### How far the match goes

| Claim | Established? | Evidence |
|---|---|---|
| our fast Null A = our line-by-line port of CellChat | **yes, exactly** | prob, p, BH q and calls identical on every TNBC row, for both versions below |
| our permutations = the ones R draws for `set.seed(1)` | **yes** | `RRandom` reproduces 5 known R outputs |
| the port = CellChat's R source | **by reading the code** | every formula and operation order follows R/modeling.R of the version named |
| our numbers = what CellChat prints when run in R | **not checked** | CellChat in R was never run (the project is Python only) |

### Two CellChat versions: the 2021 paper's (1.0.0) and the current one (2.2.0.9001)

`null_a(..., version="1.0")` and `run_native(..., version="1.0")` follow `computeCommunProb`
of sqjin/CellChat@v1.0.0, the version the Nature Communications 2021 paper describes. The
shuffles, the Hill function, agonists and antagonists, and `p = 1` where `prob = 0` are the same.
What differs:

| | 1.0.0 (paper) | 2.2.0.9001 (current, default) |
|---|---|---|
| complexes, co-receptors | per cell (geometric mean of the subunits in each cell), then trimean per group | trimean per gene, then geometric mean of the group averages |
| quantile | `stats::quantile`: `(1 - h) a + h b`, only where `b != a` | `collapse::fquantile`: `a + h (b - a)` |
| p-value | `#(prob_perm >= prob) / nboot`: ties count against | `#(prob_perm > prob) / nboot` |
| default data | `raw.use = FALSE`: data first smoothed over a protein network (`projectData`); **we run 1.0.0 on the raw data** (`raw.use = TRUE`) | raw data |

Both are run on the same CellChatDB (v2) and the same L-R pairs, so only Null A and the score
differ. The 1.0.0 port is also checked against the vectorised Null A, and matches exactly
(TNBC: 1,757,690 / 1,757,690 prob, 14,031 / 14,031 p identical).

TNBC, all 42,512 cells (`results/null_a_match_tnbc_wu2021_full.json`):

| | 2.2.0.9001 | 1.0.0 |
|---|---|---|
| rows with prob > 0 | 14,872 | 14,031 (per-cell complexes are 0 whenever one subunit is) |
| prob identical between versions | 1,753,488 of 1,757,690 rows; the 4,202 that differ involve complexes or co-receptors (1,018 of the 2,090 pairs have one) | |
| p identical, rows scored in both (14,031) | 13,765 | |
| p-value of 1.0.0 ≥ 2.2's | 13,830 of 14,031 (the ≥ rule counts ties against) | |
| calls at p < 0.05 | 12,371 in both, 653 only in 2.2, 99 only in 1.0.0 | |
| Open Problems AUPRC / odds ratio | **0.350** / 3.45 | 0.291 / 1.95 |
| Tabula Sapiens, full permutation: p < 0.05, scored rows | 15.3% (8.7% without the filter) | 15.3% (9.0% without the filter) |
| BH false calls on shuffled labels (CellChat p / (1 + #)/(B + 1)) | 835 / 0 | 772 / 0 |

So the version changes individual p-values (mostly for complexes and through ties) and changes
the ground-truth ranking noticeably, but not the main finding: in both versions, CellChat's
reported p-values are anti-conservative by about 3 times under full permutation, and its
p = 0 lets BH make false calls.

What could still differ from running CellChat in R, and why it should not:

- `log` / `exp`: numpy matches the C library bit for bit (checked on 10⁶ values each).
- R's `mean` uses long double; on Apple Silicon, long double is double, which the port assumes.
  On x86 Linux, R's `mean` may differ in the last bit. That can only flip exact ties between a
  permuted and an observed probability.
- collapse is compiled C: if the compiler fused `a + h * (b - a)` into one FMA instruction, a
  quantile could differ in the last bit. Same consequence.
- presto's p-values at exactly the 0.05 boundary can only change which L-R pairs enter, not Null A.

## Results (TNBC, Open Problems `tnbc_wu2021`, all 42,512 cells, 29 groups)

`results/null_a_match_tnbc_wu2021_full.json` (and `_quick` for 8,000 cells).

CellChat's pipeline: 1,147 CellChatDB genes in the data, 1,143 over-expressed, 2,090 L-R pairs,
29 × 29 × 2,090 = 1,757,690 (pair, source, target) rows, of which 14,872 have prob > 0. The
trimean is 0 unless more than ~25% of a group expresses the gene, so most rows are 0 and get p = 1.

| Check (14,872 testable rows) | Result |
|---|---|
| R's generator reproduces R (5 reference outputs) | all identical |
| prob, ours vs CellChat (all 1,757,690 rows) | **identical**, max difference 0 |
| p-value, same permutations (nboot = 100, seed.use = 1) | **14,872 / 14,872 identical** |
| BH q-value | **14,872 / 14,872 identical** |
| calls at p < 0.05 (CellChat's rule) | 13,024 both, 0 only ours, 0 only CellChat |
| calls at BH q < 0.05 | 13,024 both, 0 / 0 |
| independent permutations, B = 100 | Spearman 0.974, mean abs. difference in p 0.0024, 99.8% within 2 Monte-Carlo sd; calls 13,002 both, 17 only ours, 22 only CellChat |
| independent permutations, B = 1,000 | Spearman 0.948, mean abs. difference in p 0.0020, 99.6% within 2 sd (100% within 3); calls 13,007 / 18 / 17 |

Spearman is below 1 for independent draws because 12,862 of the 14,872 p-values are exactly 0
(ties), not because the null differs; the Monte-Carlo test is the relevant check.

What CellChat's own significance looks like on TNBC: 13,024 of 14,872 testable interactions
(88%) have p < 0.05, and 12,862 have p = 0 (no permutation reached the observed probability).
Because p = 0 is allowed, BH q-values are also 0 and BH removes nothing (13,024 calls either
way). CellChat reports no FDR, and its p-value cannot be read as one.

Run time on a laptop (Apple Silicon, 8 GB): the literal port 59 s, our Null A 26 s with
B = 100 and 176 s with B = 1,000.

## Evaluation against ground truth (plan step 1, TNBC cytokine truth)

`python benchmarks/run_tool.py --tool cellchat --source native --B 100 --seed 1 --nulls a`
→ `benchmarks/results/tools/cellchat_native_full.json`. Open Problems ligand-target protocol.

| | AUPRC | odds ratio (top 5%) |
|---|---|---|
| CellChat (this port), score only | **0.350** | 3.45 |
| CellChat, score kept where Null A p ≤ 0.05 | 0.350 | 3.45 |
| LIANA's CellChat re-implementation (`benchmarks/results/tnbc/results_full.json`) | 0.235 | 3.09 |
| random | 0.214 | 1.51 |

The original CellChat model (its database, cofactors and over-expression filter) ranks the
cytokine truth much better than LIANA's re-implementation. Gating by Null A changes nothing,
because the top-scoring interactions all have p < 0.05.

BH over every row (1,757,690, including prob = 0 rows with p = 1): CellChat's p-values give
12,862 calls (exactly its p = 0 rows), the `(1 + #) / (B + 1)` p-values give 0. With B = 100 the
smallest p is 1/101, above BH's threshold for that many tests. What BH "finds" here is decided by
the p = 0 convention and the number of rows counted, not by the data.

## Calibration (plan step 2, Tabula Sapiens, full permutation)

`python -m tools.cellchat.calibrate_null_a` → `results/calibration_tabula_sapiens_mammary.json`.
Tabula Sapiens - Mammary (30,934 cells, 16 cell types with ≥ 5 cells, 2 donors), raw counts →
ln(1 + CP10k). Every cell's label is shuffled once (5 shuffles), so no interaction is real, and
CellChat runs end to end with its defaults (nboot = 100, seed.use = 1).

| Fully shuffled labels | P(p < 0.05), scored rows | P(p < 0.05), all rows | BH false calls, CellChat p | BH false calls, (1 + #)/(B + 1) |
|---|---|---|---|---|
| CellChat pipeline | **15.3%** (range 10.9-17.7%) | 0.47% | 835 | 0 |
| without the over-expression filter | **8.7%** (6.5-11.1%) | 0.33% | 1,838 | 0 |
| scored rows where no permutation reached prob = 0 | **4.92%** | | | |
| expected for a calibrated test | 4.95% | ≤ 5% | 0 | 0 |

With the real labels, 7,112 of 8,952 scored rows (79%) have p < 0.05.

What this shows:
1. **Our Null A is correct.** Where nothing is selected, it is exactly calibrated (4.92% vs
   4.95%), and over all rows it is conservative, as a permutation test must be.
2. **CellChat's reported p-values are anti-conservative, about 3 times nominal**, because two
   choices made with the same labels select the rows that are reported:
   - the over-expression filter keeps genes that look high in some group under these labels
     (15.3% → 8.7% without it);
   - only rows with prob > 0 are reported, and a permutation that drives prob to 0 always loses
     to an observed prob > 0 (11.5% where that happens, 4.9% where it does not).
3. **CellChat's p = 0 convention breaks BH.** With p = 0 allowed, BH calls 835 interactions on
   data with no signal at all; with `(1 + #) / (B + 1)` it calls none.

This is the calibration baseline the other nulls are compared with. Full permutation separates
a valid test from an invalid one here only because CellChat's selection is label-dependent;
for abundance confounding, the fake-L-R-resource test is still needed (docs/PLAN.md).

## Checklist

- [x] pin the CellChat version and record it here
- [x] keep all rows (prob = 0 rows get p = 1), not only significant ones, so nulls can be compared
- [x] complexes keep CellChat's names in `ligand` / `receptor`; genes are in `ligand_complex` / `receptor_complex`
- [x] `nulls.null_a`: reproduce CellChat's own permutation test; compare `p_null_a` with `p_native`
- [ ] `nulls.null_b`

## Run the tests and the benchmark

```bash
pytest tools/cellchat
python benchmarks/run_tool.py --tool cellchat --source native --quick --B 100
```
