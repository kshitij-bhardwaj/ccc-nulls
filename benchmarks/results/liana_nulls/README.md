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
