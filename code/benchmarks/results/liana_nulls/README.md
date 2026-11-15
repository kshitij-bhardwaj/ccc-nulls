# Null A vs Null B in the LIANA+ pipeline: results

Plan: `benchmarks/liana_nulls/PLAN.md` (SHA-256 `8ba57470dcababf0…`, fixed before any result). Run: `python -m benchmarks.liana_nulls`.

## T0 fidelity and the three columns

Column 1 is the tool itself (CellChat 2.1.2 in R; CellPhoneDB 5.0.1). Column 2 is our Python Null A on the tool's own permutations. Column 3 is Null B. Calls: the tool's own rule for columns 1 and 2, BH q ≤ 0.05 over the same rows for Null A and Null B.

| Dataset | Tool | Rows where p (col. 1) = p (col. 2) | Calls, tool's rule: col. 1 / col. 2 | Null A (BH) | Null B (BH) | Null B outside Null A | Null A only: ligand / receptor fails | LIANA+ built-in vs tool rule (Jaccard) |
|---|---|---|---|---|---|---|---|---|
| tnbc_wu2021 | cellchat | 1,757,690 / 1,757,690 | 13,024 / 13,024 | 12,972 | 8,837 | 0 | 1,565 / 2,570 | 0.986 |
| tnbc_wu2021 | cellphonedb | 2,044,471 / 2,044,471 | 25,460 / 25,460 | 14,765 | 9,528 | 0 | 2,542 / 2,695 | 0.974 |
| allen_brain_atlas | cellchat | 988,701 / 988,701 | 29,458 / 29,458 | 27,973 | 13,558 | 0 | 6,622 / 7,793 | 0.972 |
| allen_brain_atlas | cellphonedb | 1,242,621 / 1,242,621 | 67,362 / 67,362 | 46,928 | 15,121 | 0 | 17,587 / 14,220 | nan |
| tabula_sapiens_mammary | cellchat | 508,160 / 508,160 | 7,112 / 7,112 | 6,873 | 3,297 | 0 | 1,299 / 2,277 | 0.948 |
| tabula_sapiens_mammary | cellphonedb | 727,552 / 727,552 | 10,637 / 10,637 | 6,201 | 2,890 | 0 | 1,648 / 1,663 | 0.984 |
| tabula_sapiens_liver | cellchat | 1,251,072 / 1,251,072 | 13,366 / 13,366 | 12,798 | 7,617 | 0 | 3,235 / 1,946 | 0.960 |
| tabula_sapiens_liver | cellphonedb | 1,653,120 / 1,653,120 | 23,558 / 23,558 | 15,773 | 7,203 | 0 | 4,032 / 4,538 | 0.991 |
| moffitt2018_hypothalamus (Amendment 1) | cellchat | 117,888 / 117,888 | 1,213 / 1,213 | 1,204 | 928 | 0 | 162 / 114 | 0.996 |
| moffitt2018_hypothalamus (Amendment 1) | cellphonedb | 156,544 / 156,544 | 3,482 / 3,482 | 2,600 | 1,350 | 0 | 951 / 299 | nan |

## T1 calibration (Tabula Sapiens, known-null data, α = 0.05)

| Tool | Tissue | Receptor genes shuffled: FDR Null A / Null B | Ligand genes shuffled: FDR Null A / Null B | Shuffled labels: P(any call) A / B |
|---|---|---|---|---|
| cellchat | mammary | 0.271 / **0.012** [0.008, 0.017] | 0.717 / **0.000** [0.000, 0.000] | 0.00 / 0.00 |
| cellchat | liver | 0.643 / **0.031** [0.019, 0.042] | 0.746 / **0.029** [0.013, 0.045] | 0.00 / 0.00 |
| cellphonedb | mammary | 1.000 / **0.000** [0.000, 0.000] | 0.970 / **0.000** [0.000, 0.000] | 0.00 / 0.00 |
| cellphonedb | liver | 0.993 / **0.000** [0.000, 0.000] | 0.955 / **0.000** [0.000, 0.000] | 0.00 / 0.00 |

## T2 power on planted edges (Tabula Sapiens Mammary, labels shuffled, 30 planted pairs)

| Tool | Fold change | Tool's rule: power / FDR | Null A: power / FDR | Null B: power / FDR | One-sided false calls A / B | Null B calls (BH) |
|---|---|---|---|---|---|---|
| cellchat | 1.5 | 0.99 / 0.879 | 0.99 / 0.844 [0.815, 0.872] | 0.19 / 0.000 [0.000, 0.000] | 251.6 / 0.0 | 9.6 |
| cellchat | 2.0 | 1.00 / 0.880 | 1.00 / 0.848 [0.819, 0.876] | 0.20 / 0.000 [0.000, 0.000] | 263.6 / 0.0 | 10.0 |
| cellchat | 3.0 | 1.00 / 0.880 | 1.00 / 0.848 [0.822, 0.875] | 0.20 / 0.000 [0.000, 0.000] | 268.0 / 0.0 | 10.0 |
| cellchat | 4.0 | 1.00 / 0.880 | 1.00 / 0.849 [0.823, 0.876] | 0.20 / 0.000 [0.000, 0.000] | 269.8 / 0.0 | 10.0 |
| cellphonedb | 1.5 | 1.00 / 0.968 | 1.00 / 0.941 [0.935, 0.948] | 0.00 / 0.000 [0.000, 0.000] | 715.4 / 0.0 | 0.0 |
| cellphonedb | 2.0 | 1.00 / 0.968 | 1.00 / 0.944 [0.937, 0.951] | 0.00 / 0.000 [0.000, 0.000] | 754.4 / 0.0 | 0.0 |
| cellphonedb | 3.0 | 1.00 / 0.968 | 1.00 / 0.945 [0.938, 0.952] | 0.00 / 0.000 [0.000, 0.000] | 768.0 / 0.0 | 0.0 |
| cellphonedb | 4.0 | 1.00 / 0.967 | 1.00 / 0.945 [0.938, 0.952] | 0.00 / 0.000 [0.000, 0.000] | 768.6 / 0.0 | 0.0 |

Replicates per cell: 5; B = 1000.

**Vacuous (Amendment 2):** Null B (BH) made no calls at cellphonedb f=1.5, cellphonedb f=2.0, cellphonedb f=3.0, cellphonedb f=4.0; its FDR there is 0 by construction (with B = 1,000 the p-value floor is 1/1001, and BH needs about m/50 rows at the floor).

### T2 under each tool's own unadjusted rule (Amendment 2)

| Tool | Fold change | Null A, tool's rule: power / FDR [95% CI] | Null B, tool's rule: power / FDR [95% CI] | One-sided false calls A / B |
|---|---|---|---|---|
| cellchat | 1.5 | 0.99 / 0.880 [0.870, 0.890] | 0.99 / 0.336 [0.293, 0.379] | 257.8 / 15.8 |
| cellchat | 2.0 | 1.00 / 0.881 [0.872, 0.889] | 1.00 / 0.320 [0.283, 0.358] | 268.4 / 14.6 |
| cellchat | 3.0 | 1.00 / 0.881 [0.874, 0.888] | 1.00 / 0.313 [0.284, 0.342] | 271.6 / 13.8 |
| cellchat | 4.0 | 1.00 / 0.881 [0.874, 0.887] | 1.00 / 0.302 [0.274, 0.330] | 273.0 / 12.6 |
| cellphonedb | 1.5 | 1.00 / 0.968 [0.964, 0.971] | 1.00 / 0.575 [0.484, 0.667] | 753.4 / 32.2 |
| cellphonedb | 2.0 | 1.00 / 0.967 [0.964, 0.971] | 1.00 / 0.558 [0.480, 0.636] | 766.0 / 28.2 |
| cellphonedb | 3.0 | 1.00 / 0.967 [0.964, 0.971] | 1.00 / 0.548 [0.455, 0.641] | 772.8 / 26.4 |
| cellphonedb | 4.0 | 1.00 / 0.967 [0.963, 0.971] | 1.00 / 0.545 [0.450, 0.641] | 773.0 / 25.8 |

## T3 biology

Primary: among Null A's calls (s ≠ t), the odds that a call Null B keeps sits on a colocalised (s, t) pair, against a call it drops; 95% CI from 2,000 bootstrap resamples of (s, t) pairs.

| Dataset | Tool | Kept: colocalised / calls | Dropped: colocalised / calls | Odds ratio [95% CI] |
|---|---|---|---|---|
| tnbc_wu2021 | cellchat | 105 / 8,345 (1.3%) | 48 / 4,009 (1.2%) | 1.05 [0.65, 1.53] |
| tnbc_wu2021 | cellphonedb | 103 / 9,052 (1.1%) | 72 / 5,067 (1.4%) | 0.80 [0.59, 1.32] |
| allen_brain_atlas | cellchat | 544 / 12,648 (4.3%) | 514 / 13,809 (3.7%) | 1.16 [0.88, 1.57] |
| allen_brain_atlas | cellphonedb | 646 / 14,113 (4.6%) | 1,264 / 30,196 (4.2%) | 1.10 [0.78, 1.48] |

Secondary (AUPRC; Null B − other, paired bootstrap 95% CI):

| Dataset | Tool | Truth | Null A | Null B | Tool's rule | Null B − Null A |
|---|---|---|---|---|---|---|
| tnbc_wu2021 | cellchat | pair_level_auprc | 0.013 | 0.013 | 0.013 | +0.000 [-0.001, 0.002] |
| tnbc_wu2021 | cellchat | openproblems_source_target | 0.013 | 0.012 | 0.013 | -0.001 [-0.006, 0.002] |
| tnbc_wu2021 | cellchat | cytokine_ligand_target | 0.350 | 0.349 | 0.350 | -0.000 [-0.032, 0.032] |
| tnbc_wu2021 | cellphonedb | pair_level_auprc | 0.013 | 0.012 | 0.013 | -0.001 [-0.010, 0.001] |
| tnbc_wu2021 | cellphonedb | openproblems_source_target | 0.010 | 0.010 | 0.010 | +0.000 [-0.003, 0.002] |
| tnbc_wu2021 | cellphonedb | cytokine_ligand_target | 0.278 | 0.242 | 0.278 | -0.036 [-0.078, 0.003] |
| allen_brain_atlas | cellchat | pair_level_auprc | 0.039 | 0.046 | 0.039 | +0.007 [0.001, 0.016] |
| allen_brain_atlas | cellchat | openproblems_source_target | 0.059 | 0.057 | 0.059 | -0.001 [-0.022, 0.014] |
| allen_brain_atlas | cellphonedb | pair_level_auprc | 0.043 | 0.049 | 0.043 | +0.006 [-0.001, 0.015] |
| allen_brain_atlas | cellphonedb | openproblems_source_target | 0.070 | 0.070 | 0.070 | +0.001 [-0.010, 0.009] |

## Decision (the plan's rule, applied as written)

**Null B is valid, but not biologically distinguishable from Null A on these truths.**

```
{
 "t1_null_b_upper_ci_le_0.05": true,
 "t2_null_b_upper_ci_le_0.05": true,
 "t2_vacuous_at": [
  "cellphonedb f=1.5",
  "cellphonedb f=2.0",
  "cellphonedb f=3.0",
  "cellphonedb f=4.0"
 ],
 "t3_primary_ci_above_1_somewhere": false,
 "t3_primary_ci_below_1_nowhere": true,
 "t3_no_secondary_significantly_worse": true,
 "worse": []
}
```

## Amendment 1 (added after the plan; `AMENDMENT_1.md`, SHA-256 `85e3e8f7811807cb…`)

Amendment 2 (`AMENDMENT_2.md`, SHA-256 `bcfff59ab5eb1ccc…`) adds the T2 tool-rule table above. Not part of the pre-registered decision above. The amendment was written before any T3 result of this pipeline; it was not changed after a smoke test of A2 on the Allen brain.

### A1. Hypothalamus (Moffitt 2018): MERFISH neighbourhood truth, 10 of 72 ordered class pairs colocalised

| Tool | Kept: colocalised / calls | Dropped: colocalised / calls | Odds ratio [95% CI] |
|---|---|---|---|
| cellchat | 164 / 696 (23.6%) | 33 / 235 (14.0%) | 1.89 [1.18, 3.34] |
| cellphonedb | 147 / 1,042 (14.1%) | 210 / 1,025 (20.5%) | 0.64 [0.35, 1.60] |

### A2. Downstream response of the receiver (Null A's calls, s ≠ t; responds = z ≥ 1.645)

D1: odds of response, calls Null B keeps vs drops. D2: kept vs dropped because the receptor side fails (predicted > 1) or the ligand side fails (no prediction). Bootstrap over (footprint, receiver).

| Dataset | Tool | Truth | Calls in scope | Responding: kept / dropped | D1 | D1 (z ≥ 1) | D2 receptor | D2 ligand |
|---|---|---|---|---|---|---|---|---|
| tnbc_wu2021 | cellchat | cytosig | 212 | 8.1% / 18.5% | 0.39 [0.05, 8.66] | 3.58 [0.51, 72.85] | 0.39 [0.05, 8.66] | – |
| tnbc_wu2021 | cellchat | progeny | 531 | 17.3% / 1.8% | 11.69 [2.90, 63.73] | 12.46 [4.56, 68.31] | 11.69 [2.90, 63.73] | – |
| tnbc_wu2021 | cellphonedb | cytosig | 2,043 | 8.9% / 4.3% | 2.16 [1.27, 3.71] | 1.46 [1.08, 2.18] | 4.12 [1.11, 67.23] | 1.48 [0.93, 3.18] |
| tnbc_wu2021 | cellphonedb | progeny | 413 | 22.3% / 9.7% | 2.69 [1.73, 4.28] | 3.15 [2.08, 5.29] | 81.99 [20.78, 187.37] | 0.66 [0.45, 0.84] |
| allen_brain_atlas | cellchat | cytosig | 906 | 8.7% / 1.7% | 5.38 [1.86, 12.10] | 1.68 [0.84, 5.20] | 9.67 [1.84, 55.69] | 1.95 [0.62, 9.70] |
| allen_brain_atlas | cellchat | progeny | 1,145 | 5.2% / 5.8% | 0.88 [0.06, 2.39] | 1.02 [0.54, 1.70] | 0.90 [0.03, 27.40] | 0.85 [0.49, 5.99] |
| allen_brain_atlas | cellphonedb | cytosig | 1,418 | 11.6% / 4.9% | 2.55 [1.55, 3.89] | 1.76 [1.18, 2.84] | 6.68 [1.84, 39.17] | 1.06 [0.65, 1.66] |
| allen_brain_atlas | cellphonedb | progeny | 308 | 5.4% / 8.1% | 0.65 [0.02, 17.43] | 0.57 [0.03, 1.74] | 0.38 [0.01, 10.44] | 7.32 [0.36, 19.24] |
| tabula_sapiens_mammary | cellchat | cytosig | 395 | 6.9% / 1.0% | 7.04 [1.16, 51.16] | 1.63 [0.79, 3.72] | 4.74 [0.77, 36.88] | 9.56 [2.47, 19.59] |
| tabula_sapiens_mammary | cellchat | progeny | 419 | 26.3% / 4.9% | 6.94 [3.51, 16.53] | 4.92 [2.04, 12.68] | 6.76 [3.40, 16.20] | 5.41 [0.65, 16.23] |
| tabula_sapiens_mammary | cellphonedb | cytosig | 1,299 | 5.9% / 5.3% | 1.12 [0.38, 3.50] | 1.60 [1.04, 2.66] | 1.26 [0.30, 13.36] | 0.88 [0.36, 2.39] |
| tabula_sapiens_mammary | cellphonedb | progeny | 213 | 26.9% / 17.8% | 1.70 [1.08, 2.94] | 1.86 [1.17, 3.35] | 9.55 [1.90, 57.65] | 1.04 [0.66, 1.52] |
| tabula_sapiens_liver | cellchat | cytosig | 813 | 10.7% / 15.5% | 0.65 [0.23, 2.69] | 1.13 [0.56, 2.48] | 0.59 [0.17, 5.91] | 0.94 [0.47, 4.04] |
| tabula_sapiens_liver | cellchat | progeny | 399 | 10.2% / 4.5% | 2.39 [0.68, 10.68] | 1.30 [0.58, 3.07] | 2.20 [0.42, 26.60] | 2.77 [1.39, 7.74] |
| tabula_sapiens_liver | cellphonedb | cytosig | 2,483 | 9.4% / 8.8% | 1.08 [0.57, 2.40] | 1.10 [0.76, 1.63] | 0.99 [0.48, 2.70] | 1.46 [0.85, 4.45] |
| tabula_sapiens_liver | cellphonedb | progeny | 319 | 9.5% / 10.6% | 0.89 [0.69, 1.25] | 1.23 [0.74, 2.06] | 9.41 [1.74, 27.27] | 0.61 [0.50, 0.74] |
| moffitt2018_hypothalamus | cellchat | cytosig | 5 | – | – | – | – | – |
| moffitt2018_hypothalamus | cellchat | progeny | 59 | 3.8% / 0.0% | 0.74 [0.06, 2.71] | 0.32 [0.04, 2.42] | 0.25 [0.01, 1.00] | 0.54 [0.03, 2.14] |
| moffitt2018_hypothalamus | cellphonedb | cytosig | 62 | 9.3% / 12.5% | 0.71 [0.09, 3.05] | 1.38 [0.26, 6.28] | 0.71 [0.09, 3.05] | – |
| moffitt2018_hypothalamus | cellphonedb | progeny | 5 | – | – | – | – | – |

### Extended verdict (Amendment 1)

- Criterion (2) over all three spatial truths: CI above 1 somewhere = True, below 1 somewhere = False → **met**.
- Downstream truths: **against**, from: tnbc_wu2021/cellphonedb/progeny/D2_ligand 0.66 [0.45, 0.84]; tabula_sapiens_liver/cellphonedb/progeny/D2_ligand 0.61 [0.50, 0.74].

Note: the T2 part of criterion (1) holds only vacuously at cellphonedb f=1.5, cellphonedb f=2.0, cellphonedb f=3.0, cellphonedb f=4.0 (Amendment 2).

## Amendment 3 (an addition; not part of the pre-registered plan)

`AMENDMENT_3.md` (SHA-256 `21fe25fc0127b0d2…`), written after every earlier result and before any of these. It changes neither verdict above.

### A3.1 T2 at small effects and B = 5,000

Check: f = 1.5, B = 1,000 reproduces `power.json` (Null A and Null B: calls, power, FDP, one-sided calls) in 10 of 10 runs.

Power / FDR, mean over replicates (FDR 95% t-CI for Null B).

| Tool | Setting | Null A (BH) | Null B (BH) | Null B calls (BH) | Null B, tool's rule | B2 (BH) |
|---|---|---|---|---|---|---|
| cellchat | f=1.1, B=1000 | 0.00 / 0.00 | 0.00 / 0.00 [0.00, 0.00] | 0.0 | 0.03 / 0.90 | 0.00 / 0.00 |
| cellchat | f=1.2, B=1000 | 0.20 / 0.91 | 0.00 / 0.00 [0.00, 0.00] | 0.0 | 0.19 / 0.69 | 0.05 / 0.00 |
| cellchat | f=1.2, B=5000 | 0.19 / 0.91 | 0.00 / 0.00 [0.00, 0.00] | 0.0 | 0.19 / 0.69 | 0.05 / 0.00 |
| cellchat | f=1.3, B=1000 | 0.76 / 0.85 | 0.00 / 0.00 [0.00, 0.00] | 0.0 | 0.76 / 0.39 | 0.33 / 0.00 |
| cellchat | f=1.5, B=1000 | 0.99 / 0.84 | 0.19 / 0.00 [0.00, 0.00] | 9.6 | 0.99 / 0.34 | 0.49 / 0.00 |
| cellchat | f=1.5, B=5000 | 0.99 / 0.84 | 0.97 / 0.02 [-0.03, 0.06] | 50.0 | 0.99 / 0.33 | 0.49 / 0.00 |
| cellphonedb | f=1.1, B=1000 | 0.27 / 0.34 | 0.00 / 0.00 [0.00, 0.00] | 0.0 | 0.71 / 0.67 | 0.03 / 0.00 |
| cellphonedb | f=1.2, B=1000 | 0.93 / 0.92 | 0.00 / 0.00 [0.00, 0.00] | 0.0 | 0.93 / 0.62 | 0.14 / 0.00 |
| cellphonedb | f=1.2, B=5000 | 0.95 / 0.92 | 0.00 / 0.00 [0.00, 0.00] | 0.0 | 0.93 / 0.61 | 0.14 / 0.00 |
| cellphonedb | f=1.3, B=1000 | 0.99 / 0.93 | 0.00 / 0.00 [0.00, 0.00] | 0.0 | 0.99 / 0.60 | 0.30 / 0.00 |
| cellphonedb | f=1.5, B=1000 | 1.00 / 0.94 | 0.00 / 0.00 [0.00, 0.00] | 0.0 | 1.00 / 0.58 | 0.47 / 0.00 |
| cellphonedb | f=1.5, B=5000 | 1.00 / 0.94 | 0.20 / 0.00 [0.00, 0.00] | 10.4 | 1.00 / 0.58 | 0.47 / 0.00 |

**Ablation on T2** (f = 1.5, B = 1,000; power / FDR / one-sided false calls, mean of 5 replicates):

| Variant | CellChat | CellPhoneDB |
|---|---|---|
| Null A (BH) | 0.99 / 0.84 / 252 | 1.00 / 0.94 / 715 |
| **Null B (BH), full** | 0.19 / 0.00 / 0 | 0.00 / 0.00 / 0 |
| Null B, ligand side only | 0.99 / 0.73 / 133 | 0.97 / 0.89 / 376 |
| Null B, receptor side only | 0.99 / 0.75 / 143 | 1.00 / 0.90 / 388 |
| Null B, BY | 0.00 / 0.00 / 0 | 0.00 / 0.00 / 0 |
| Null B, tool's rule | 0.99 / 0.34 / 16 | 1.00 / 0.58 / 32 |
| Null A, tool's rule | 0.99 / 0.88 / 258 | 1.00 / 0.97 / 753 |
| Null B, testable family | – | 0.00 / 0.00 / 0 |
| Null B (#9), gate | – | 0.00 / 0.00 / 0 |
| Null B (#9), testable | – | 0.00 / 0.00 / 0 |
| Null B (#9), BY | – | 0.00 / 0.00 / 0 |
| **B2 (CAMERA), full** | 0.49 / 0.00 / 0 | 0.47 / 0.00 / 0 |
| B2, no correction (old premise) | 0.53 / 0.00 / 0 | 0.51 / 0.00 / 0 |
| B2, unmatched decoys | 0.95 / 0.00 / 0 | 0.91 / 0.00 / 0 |
| B2, joint (sides combined) | 0.59 / 0.91 / 282 | 0.59 / 0.95 / 457 |

### A3.2 CellPhoneDB's Null B as pushed in #9

#9's own code: float64 means, `default_rng(0)` permutations; ours: CellPhoneDB's permutation sequence and its float32 means.

| Dataset | Calls, ours (BH, gate) | Calls, #9 (BH, gate) | Same call | Jaccard | p identical (family) | Spearman p | #9 testable | #9 BY |
|---|---|---|---|---|---|---|---|---|
| tnbc_wu2021 | 9,528 | 9,518 | 0.9999 | 0.985 | 0.643 | 0.9919 | 0 | 8,218 |
| allen_brain_atlas | 15,121 | 15,237 | 0.9994 | 0.970 | 0.723 | 0.9815 | 0 | 0 |
| tabula_sapiens_mammary | 2,890 | 2,902 | 0.9999 | 0.983 | 0.642 | 0.9927 | 0 | 2,409 |
| tabula_sapiens_liver | 7,203 | 7,170 | 0.9999 | 0.977 | 0.640 | 0.9940 | 0 | 0 |
| moffitt2018_hypothalamus | 1,350 | 1,350 | 1.0000 | 1.000 | 0.878 | 0.9975 | 1,302 | 1,302 |

### A3.3 Null B2: decoy genes with the CAMERA correction

Validity (Tabula Sapiens, 10 replicates). Shuffled labels: side-level type-I error at α = 0.05 (95% t-CI) and P(any BH call). Receptor genes shuffled: FDR (95% t-CI). Calibrated = type-I CI reaches ≤ 0.05 and FDR upper bound ≤ 0.05.

| Tissue | Tool | Variant | Type-I (labels shuffled) | P(any call) | FDR (receptor shuffled) | Calibrated |
|---|---|---|---|---|---|---|
| tabula_sapiens_mammary | cellchat | **B2 (CAMERA), full** | 0.033 [0.030, 0.035] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_mammary | cellchat | B2, no correction (old premise) | 0.042 [0.039, 0.045] | 1/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_mammary | cellchat | B2, unmatched decoys | 0.032 [0.029, 0.036] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_mammary | cellchat | B2, joint (sides combined) | 0.034 [0.031, 0.038] | 2/10 | 0.399 [0.394, 0.404] | – |
| tabula_sapiens_mammary | cellphonedb | **B2 (CAMERA), full** | 0.032 [0.029, 0.035] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_mammary | cellphonedb | B2, no correction (old premise) | 0.041 [0.038, 0.045] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_mammary | cellphonedb | B2, unmatched decoys | 0.032 [0.029, 0.035] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_mammary | cellphonedb | B2, joint (sides combined) | 0.041 [0.034, 0.047] | 0/10 | 1.000 [1.000, 1.000] | – |
| tabula_sapiens_liver | cellchat | **B2 (CAMERA), full** | 0.033 [0.031, 0.035] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_liver | cellchat | B2, no correction (old premise) | 0.046 [0.044, 0.049] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_liver | cellchat | B2, unmatched decoys | 0.035 [0.032, 0.037] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_liver | cellchat | B2, joint (sides combined) | 0.035 [0.034, 0.036] | 0/10 | 0.654 [0.652, 0.657] | – |
| tabula_sapiens_liver | cellphonedb | **B2 (CAMERA), full** | 0.032 [0.029, 0.035] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_liver | cellphonedb | B2, no correction (old premise) | 0.045 [0.042, 0.047] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_liver | cellphonedb | B2, unmatched decoys | 0.035 [0.032, 0.038] | 0/10 | 0.000 [0.000, 0.000] | yes |
| tabula_sapiens_liver | cellphonedb | B2, joint (sides combined) | 0.035 [0.031, 0.039] | 0/10 | 1.000 [1.000, 1.000] | – |

Real data: BH calls over the tool's family (B2's gate: every subunit in ≥ 10% of its type's cells).

| Dataset | Tool | Family | In family, passing B2's gate | B2 | No correction | Unmatched | Joint | Jaccard B2 vs Null B |
|---|---|---|---|---|---|---|---|---|
| tnbc_wu2021 | cellchat | 27,942 | 19,521 | 204 | 214 | 185 | 6,256 | 0.023 |
| tnbc_wu2021 | cellphonedb | 18,815 | 13,442 | 13 | 11 | 12 | 1,609 | 0.001 |
| allen_brain_atlas | cellchat | 162,047 | 80,322 | 10 | 10 | 4 | 1,609 | 0.001 |
| allen_brain_atlas | cellphonedb | 135,279 | 37,140 | 1 | 1 | 1 | 399 | 0.000 |
| tabula_sapiens_mammary | cellchat | 71,456 | 20,237 | 66 | 71 | 71 | 2,240 | 0.020 |
| tabula_sapiens_mammary | cellphonedb | 11,767 | 6,977 | 17 | 18 | 2 | 1,032 | 0.006 |
| tabula_sapiens_liver | cellchat | 80,794 | 31,063 | 3 | 5 | 43 | 3,938 | 0.000 |
| tabula_sapiens_liver | cellphonedb | 31,315 | 18,529 | 3 | 3 | 29 | 1,843 | 0.000 |
| moffitt2018_hypothalamus | cellchat | 5,273 | 2,888 | 5 | 5 | 8 | 209 | 0.005 |
| moffitt2018_hypothalamus | cellphonedb | 5,491 | 2,011 | 0 | 0 | 0 | 77 | 0.000 |

### A3.4 Ablation on the five data sets

Calls (Jaccard with the full Null B). One part changed at a time.

**cellchat**

| Variant | tnbc_wu2021 | allen_brain_atlas | tabula_sapiens_mammary | tabula_sapiens_liver | moffitt2018_hypothalamus |
|---|---|---|---|---|---|
| Null A (BH) | 12,972 (0.68) | 27,973 (0.48) | 6,873 (0.48) | 12,798 (0.60) | 1,204 (0.77) |
| **Null B (BH), full** | 8,837 (1.00) | 13,558 (1.00) | 3,297 (1.00) | 7,617 (1.00) | 928 (1.00) |
| Null B, ligand side only | 12,069 (0.73) | 26,902 (0.50) | 6,137 (0.54) | 10,157 (0.75) | 1,119 (0.83) |
| Null B, receptor side only | 10,792 (0.82) | 25,785 (0.53) | 5,010 (0.66) | 11,908 (0.64) | 1,197 (0.78) |
| Null B, BY | 7,905 (0.89) | 0 (0.00) | 0 (0.00) | 0 (0.00) | 0 (0.00) |
| Null B, tool's rule | 8,961 (0.99) | 15,876 (0.85) | 3,820 (0.86) | 8,666 (0.88) | 939 (0.99) |
| Null A, tool's rule | 13,033 (0.68) | 29,573 (0.46) | 7,122 (0.46) | 13,375 (0.57) | 1,211 (0.77) |
| **B2 (CAMERA), full** | 204 (0.02) | 10 (0.00) | 66 (0.02) | 3 (0.00) | 5 (0.01) |
| B2, no correction (old premise) | 214 (0.02) | 10 (0.00) | 71 (0.02) | 5 (0.00) | 5 (0.01) |
| B2, unmatched decoys | 185 (0.02) | 4 (0.00) | 71 (0.02) | 43 (0.01) | 8 (0.01) |
| B2, joint (sides combined) | 6,256 (0.50) | 1,609 (0.07) | 2,240 (0.36) | 3,938 (0.34) | 209 (0.14) |

**cellphonedb**

| Variant | tnbc_wu2021 | allen_brain_atlas | tabula_sapiens_mammary | tabula_sapiens_liver | moffitt2018_hypothalamus |
|---|---|---|---|---|---|
| Null A (BH) | 14,765 (0.65) | 46,928 (0.32) | 6,201 (0.47) | 15,773 (0.46) | 2,600 (0.52) |
| **Null B (BH), full** | 9,528 (1.00) | 15,121 (1.00) | 2,890 (1.00) | 7,203 (1.00) | 1,350 (1.00) |
| Null B, ligand side only | 13,495 (0.71) | 44,199 (0.34) | 6,031 (0.48) | 15,065 (0.48) | 2,177 (0.62) |
| Null B, receptor side only | 13,982 (0.68) | 49,588 (0.30) | 6,187 (0.47) | 15,859 (0.45) | 3,384 (0.40) |
| Null B, BY | 8,243 (0.87) | 0 (0.00) | 0 (0.00) | 0 (0.00) | 1,303 (0.97) |
| Null B, tool's rule | 9,959 (0.96) | 19,181 (0.79) | 3,213 (0.90) | 8,272 (0.87) | 1,460 (0.92) |
| Null A, tool's rule | 14,870 (0.64) | 49,561 (0.31) | 6,438 (0.45) | 16,364 (0.44) | 2,631 (0.51) |
| Null B, testable family | 0 (0.00) | 13,723 (0.91) | 0 (0.00) | 0 (0.00) | 1,303 (0.97) |
| Null B (#9), gate | 9,518 (0.98) | 15,237 (0.97) | 2,902 (0.98) | 7,170 (0.98) | 1,350 (1.00) |
| Null B (#9), testable | 0 (0.00) | 0 (0.00) | 0 (0.00) | 0 (0.00) | 1,302 (0.96) |
| Null B (#9), BY | 8,218 (0.86) | 0 (0.00) | 2,409 (0.83) | 0 (0.00) | 1,302 (0.96) |
| **B2 (CAMERA), full** | 13 (0.00) | 1 (0.00) | 17 (0.01) | 3 (0.00) | 0 (0.00) |
| B2, no correction (old premise) | 11 (0.00) | 1 (0.00) | 18 (0.01) | 3 (0.00) | 0 (0.00) |
| B2, unmatched decoys | 12 (0.00) | 1 (0.00) | 2 (0.00) | 29 (0.00) | 0 (0.00) |
| B2, joint (sides combined) | 1,609 (0.14) | 399 (0.02) | 1,032 (0.20) | 1,843 (0.18) | 77 (0.04) |

Biology of each variant: odds ratio, Null A's calls the variant keeps vs drops (95% bootstrap CI; ↑ above 1, ↓ below 1). – where the variant keeps none or drops none; (k = n) where it keeps fewer than 20.

**cellchat**

| Variant | tnbc_wu2021 spatial | allen_brain_atlas spatial | moffitt2018_hypothalamus spatial | tnbc_wu2021 cytosig | tnbc_wu2021 progeny | allen_brain_atlas cytosig | allen_brain_atlas progeny | tabula_sapiens_mammary cytosig | tabula_sapiens_mammary progeny | tabula_sapiens_liver cytosig | tabula_sapiens_liver progeny | moffitt2018_hypothalamus cytosig | moffitt2018_hypothalamus progeny |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Null B (BH), full** | 1.05 [0.65, 1.53] | 1.16 [0.88, 1.57] | 1.89 [1.18, 3.34]↑ | 0.39 [0.05, 8.66] | 11.69 [2.90, 63.73]↑ | 5.38 [1.86, 12.10]↑ | 0.88 [0.06, 2.39] | 7.04 [1.16, 51.16]↑ | 6.94 [3.51, 16.53]↑ | 0.65 [0.23, 2.69] | 2.39 [0.68, 10.68] | – | 0.74 [0.06, 2.71] |
| Null B, ligand side only | 1.19 [0.68, 3.70] | 1.39 [0.93, 2.10] | 2.03 [1.13, 5.07]↑ | – | – | 0.84 [0.24, 3.95] | 0.85 [0.37, 8.12] | 6.37 [1.63, 13.08]↑ | 2.04 [0.22, 5.83] | 1.15 [0.43, 5.41] | 2.40 [1.13, 7.18]↑ | – | 0.52 [0.03, 2.05] |
| Null B, receptor side only | 0.95 [0.41, 2.45] | 0.98 [0.68, 1.49] | 1.46 [0.72, 3.63] | 0.39 [0.05, 8.66] | 11.69 [2.90, 63.73]↑ | 8.31 [1.62, 49.19]↑ | 0.90 [0.01, 27.29] | 3.24 [0.50, 26.63] | 6.31 [3.18, 14.66]↑ | 0.69 [0.19, 34.69] | 2.03 [0.38, 24.97] | – | 0.23 [0.01, 0.88]↓ |
| Null B, BY | 0.72 [0.52, 1.12] | – | – | 0.82 [0.16, 12.01] | 1.75 [0.65, 5.81] | – | – | – | – | – | – | – | – |
| Null B, tool's rule | 1.22 [0.64, 2.04] | 1.12 [0.86, 1.54] | 1.77 [1.11, 3.11]↑ | 0.37 [0.04, 8.41] | 11.69 [2.90, 63.73]↑ | 6.47 [1.56, 29.51]↑ | 1.48 [0.38, 8.65] | 3.09 [0.52, 22.47] | 4.67 [2.44, 10.15]↑ | 0.47 [0.16, 2.17] | 1.32 [0.37, 6.12] | – | 0.74 [0.06, 2.71] |
| Null A, tool's rule | – | – | – | – | – | – | – | – | – | – | – | – | – |
| **B2 (CAMERA), full** | 0.22 [0.12, 0.66]↓ | 3.43 [1.09, 15.85]↑ (k=8) | 0.74 [0.26, 5.13] (k=2) | – | 2.47 [0.63, 12.44] (k=18) | – | – | 2.06 [0.80, 10.44] (k=5) | 2.49 [0.70, 25.46] (k=1) | – | – | 3.00 [0.71, 11.00] (k=1) | – |
| B2, no correction (old premise) | 0.21 [0.11, 0.62]↓ | 3.43 [1.09, 15.85]↑ (k=8) | 0.74 [0.26, 5.13] (k=2) | – | 1.88 [0.50, 5.86] | – | – | 2.06 [0.80, 10.44] (k=5) | 2.49 [0.70, 25.46] (k=1) | 1.38 [0.42, 13.82] (k=2) | – | 3.00 [0.71, 11.00] (k=1) | – |
| B2, unmatched decoys | 0.25 [0.13, 0.75]↓ | 3.43 [1.51, 23.41]↑ (k=3) | 1.87 [0.36, 18.33] (k=3) | 0.53 [0.20, 3.56] (k=8) | 0.47 [0.12, 3.51] (k=14) | – | – | – | – | – | – | – | – |
| B2, joint (sides combined) | 0.70 [0.54, 1.02] | 1.19 [0.88, 1.61] | 1.92 [0.92, 4.36] | 1.17 [0.51, 4.71] | 1.62 [0.65, 4.59] | 1.77 [0.36, 7.02] | 0.22 [0.10, 0.97]↓ | 2.43 [0.50, 13.40] | 0.87 [0.22, 2.48] | 0.47 [0.16, 1.27] | 1.28 [0.43, 5.41] | 0.33 [0.09, 1.40] (k=4) | 1.16 [0.31, 13.01] (k=8) |

**cellphonedb**

| Variant | tnbc_wu2021 spatial | allen_brain_atlas spatial | moffitt2018_hypothalamus spatial | tnbc_wu2021 cytosig | tnbc_wu2021 progeny | allen_brain_atlas cytosig | allen_brain_atlas progeny | tabula_sapiens_mammary cytosig | tabula_sapiens_mammary progeny | tabula_sapiens_liver cytosig | tabula_sapiens_liver progeny | moffitt2018_hypothalamus cytosig | moffitt2018_hypothalamus progeny |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **Null B (BH), full** | 0.80 [0.59, 1.32] | 1.10 [0.78, 1.48] | 0.64 [0.35, 1.60] | 2.16 [1.27, 3.71]↑ | 2.69 [1.73, 4.28]↑ | 2.55 [1.55, 3.89]↑ | 0.65 [0.02, 17.43] | 1.12 [0.38, 3.50] | 1.70 [1.08, 2.94]↑ | 1.08 [0.57, 2.40] | 0.89 [0.69, 1.25] | 0.71 [0.09, 3.05] | – |
| Null B, ligand side only | 0.80 [0.47, 2.11] | 1.04 [0.74, 1.47] | 0.74 [0.33, 2.70] | 1.39 [0.88, 3.39] | 0.33 [0.19, 0.47]↓ | 0.50 [0.28, 0.98]↓ | 10.96 [2.45, 23.14]↑ | 0.76 [0.23, 3.11] | 0.55 [0.31, 0.88]↓ | 1.54 [0.72, 7.90] | 0.93 [0.69, 1.16] | – | – |
| Null B, receptor side only | 1.04 [0.66, 1.66] | 1.12 [0.82, 1.65] | 0.69 [0.32, 1.74] | 3.52 [0.93, 58.36] | 91.14 [22.96, 209.01]↑ | 9.62 [3.47, 37.81]↑ | 0.25 [0.01, 6.95] | 1.30 [0.32, 13.01] | 9.33 [1.88, 56.25]↑ | 0.88 [0.40, 2.36] | 8.56 [1.41, 25.10]↑ | 0.71 [0.09, 3.05] | – |
| Null B, BY | 0.60 [0.34, 1.17] | – | 0.66 [0.35, 1.75] | 3.20 [1.88, 5.24]↑ | 1.77 [1.24, 2.81]↑ | – | – | – | – | – | – | 0.71 [0.09, 3.05] | – |
| Null B, tool's rule | 0.94 [0.70, 1.49] | 1.26 [0.98, 1.63] | 0.75 [0.50, 1.52] | 2.19 [1.26, 3.67]↑ | 2.54 [1.62, 4.04]↑ | 2.93 [2.10, 4.42]↑ | 0.47 [0.02, 13.18] | 1.10 [0.39, 3.86] | 1.65 [1.01, 2.94]↑ | 0.95 [0.49, 2.26] | 1.63 [1.18, 2.26]↑ | 0.71 [0.09, 3.05] | – |
| Null A, tool's rule | – | – | – | – | – | – | – | – | – | – | – | – | – |
| Null B, testable family | – | 1.11 [0.77, 1.52] | 0.66 [0.35, 1.75] | – | – | 2.91 [1.73, 4.50]↑ | 0.53 [0.03, 3.88] | – | – | – | – | 0.71 [0.09, 3.05] | – |
| Null B (#9), gate | 0.78 [0.59, 1.32] | 1.10 [0.78, 1.50] | 0.64 [0.35, 1.60] | 2.11 [1.24, 3.63]↑ | 2.69 [1.73, 4.28]↑ | 2.94 [1.78, 4.55]↑ | 0.65 [0.02, 17.43] | 1.06 [0.37, 3.35] | 1.70 [1.08, 2.94]↑ | 1.08 [0.57, 2.40] | 1.21 [0.91, 1.72] | 0.71 [0.09, 3.05] | – |
| Null B (#9), testable | – | – | 0.67 [0.36, 1.75] | – | – | – | – | – | – | – | – | 1.22 [0.18, 4.97] | – |
| Null B (#9), BY | 0.61 [0.31, 1.36] | – | 0.67 [0.36, 1.75] | 3.04 [1.79, 4.97]↑ | 3.13 [2.01, 5.50]↑ | – | – | 1.22 [0.41, 3.91] | 1.70 [1.08, 2.94]↑ | – | – | 1.22 [0.18, 4.97] | – |
| **B2 (CAMERA), full** | 2.94 [1.39, 10.35]↑ (k=13) | 66.61 [15.72, 206.71]↑ (k=1) | – | – | 50.04 [31.84, 69.35]↑ (k=4) | 4.47 [1.60, 22.87]↑ (k=1) | – | 1.90 [0.77, 11.80] (k=4) | – | 3.35 [1.10, 17.62]↑ (k=1) | – | – | – |
| B2, no correction (old premise) | 3.45 [1.63, 12.39]↑ (k=11) | 66.61 [15.72, 206.71]↑ (k=1) | – | – | – | 4.47 [1.60, 22.87]↑ (k=1) | – | 1.56 [0.65, 7.40] (k=5) | – | 3.35 [1.10, 17.62]↑ (k=1) | – | – | – |
| B2, unmatched decoys | 3.18 [1.57, 10.73]↑ (k=12) | 7.40 [2.49, 34.64]↑ (k=1) | – | 4.35 [1.44, 27.16]↑ (k=1) | – | – | – | – | – | – | – | – | – |
| B2, joint (sides combined) | 0.82 [0.39, 1.27] | 1.69 [0.94, 2.93] | 1.49 [0.68, 3.41] | 2.05 [0.35, 6.75] | 39.49 [6.90, 358.85]↑ | 0.28 [0.08, 1.35] | 0.75 [0.22, 4.39] (k=19) | 1.02 [0.20, 3.49] | 0.45 [0.24, 0.74]↓ | 0.55 [0.12, 1.30] | 25.93 [12.50, 73.00]↑ | 0.34 [0.10, 6.74] (k=10) | – |

Notes on reading Amendment 3:

- "B2, no correction" is B2 without the CAMERA variance term. It is not the October decoy null, which used a different statistic and was not rerun here.
- Variants with 0 calls under BH or BY over a large family (#9's testable family, BY on Allen and Tabula Sapiens) are the permutation floor: the smallest p is 1/(B+1), and too few rows reach it to pass the step-up threshold.
- Two display rules were set after seeing the results; they change no number. An odds ratio is "–" when the variant keeps none or drops none of Null A's calls (the 0.5 correction otherwise prints a number), and carries its kept count when that count is below 20.

## Status

| Component | Plan | Result files | Status |
|---|---|---|---|
| T0 fidelity + columns 1-3 | pre-registered | 4/4 | done |
| T1 calibration | pre-registered | 1/1 | done |
| T2 power | pre-registered | 1/1 | done |
| T2 under the tools' own rule | Amendment 2 | 1/1 | done |
| T3 biology (TNBC, Allen brain) | pre-registered | 1/1 | done |
| A1 hypothalamus: core | Amendment 1 | 1/1 | done |
| A1 hypothalamus: spatial truth | Amendment 1 | 1/1 | done |
| A2 downstream truths (CytoSig, PROGENy) | Amendment 1 | 1/1 | done |
| A3.1 T2 small effects, B = 5,000 | Amendment 3 | 1/1 | done |
| A3.2 #9's CellPhoneDB Null B | Amendment 3 | 5/5 | done |
| A3.3 B2 on real data | Amendment 3 | 5/5 | done |
| A3.3 B2 validity | Amendment 3 | 1/1 | done |
| A3.4 ablation | Amendment 3 | 1/1 | done |
