# CellChat

One owner, assigned within the team. Only the owner edits this folder; anything shared goes in
`src/ccc_nulls/` through a reviewed pull request.

| | |
|---|---|
| Tool | CellChat **2.2.0.9001**, `jinworks/CellChat@75253cd` (latest on GitHub, March 2026), **ported to Python**: no R needed |
| Statistic tested by the nulls | `prob`, CellChat's communication probability (LIANA's re-implementation calls it `lr_probs`) |
| Native null (Null A) | Yes: `computeCommunProb` permutes the cell labels, `nboot = 100`, `seed.use = 1` |
| Extra Python dependency | `rdata` (reads CellChatDB from the package's `.rda` once): `pip install rdata` |

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
