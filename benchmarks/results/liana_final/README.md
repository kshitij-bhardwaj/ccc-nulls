# Our nulls inside LIANA+'s consensus: results

Plan: `benchmarks/liana_final/PLAN.md` (SHA-256 `b600baf70d17313a…`), amended by `AMENDMENT_1.md` (SHA-256 `13396316c60a5f92…`): F1 and F2 judged on the ranking, after one power replicate showed the call rule could not separate the variants. Run: `python -m benchmarks.liana_final`.
Baseline: LIANA+'s five methods plus CellChat (liana6). A call (F0, F3) is specificity_rank ≤ 0.05.

## F0 fidelity and calls

| Dataset | Rows | liana5 = rank_aggregate | LIANA+ as shipped (5 methods) | **+ CellChat (baseline)** | **swap: Null B** | swap: B2 | 7th/8th: Null A | **7th/8th: Null B** | 7th/8th: B2 |
|---|---|---|---|---|---|---|---|---|---|
| tnbc_wu2021 | 154,213 | yes | 34,013 (0.54) | 23,466 (1.00) | 24,005 (0.69) | 34,458 (0.59) | 27,505 (0.75) | 23,122 (0.74) | 33,275 (0.60) |
| allen_brain_atlas | 235,366 | yes | 50,770 (0.65) | 42,879 (1.00) | 39,249 (0.67) | 50,911 (0.62) | 53,451 (0.74) | 42,186 (0.68) | 52,213 (0.63) |
| tabula_sapiens_mammary | 86,001 | yes | 17,153 (0.70) | 15,946 (1.00) | 13,208 (0.72) | 19,373 (0.64) | 19,093 (0.82) | 15,041 (0.73) | 18,852 (0.66) |
| tabula_sapiens_liver | 196,381 | yes | 42,982 (0.60) | 33,349 (1.00) | 33,740 (0.75) | 45,527 (0.65) | 39,560 (0.77) | 34,632 (0.75) | 44,768 (0.65) |
| moffitt2018_hypothalamus | 10,197 | yes | 766 (0.40) | 1,456 (1.00) | 1,292 (0.89) | 1,937 (0.65) | 1,757 (0.75) | 1,480 (0.69) | 1,814 (0.70) |

Calls (Jaccard with the baseline's calls).

## F1 calibration (Tabula Sapiens, known nulls; AMENDMENT_1.md)

Share of rows with a shuffled side in the top 1% of each ranking, mean of 3 replicates [95% t-CI]. The first line is the share among all rows. Lower is better.

| Variant | mammary, receptors shuffled | mammary, ligands shuffled | liver, receptors shuffled | liver, ligands shuffled |
|---|---|---|---|---|
| all rows | 0.920 | 0.974 | 0.929 | 0.963 |
| LIANA+ as shipped (5 methods) | 0.797 [0.784, 0.811] | 0.811 [0.807, 0.816] | 0.881 [0.876, 0.887] | 0.803 [0.784, 0.821] |
| **+ CellChat (baseline)** | 0.799 [0.782, 0.815] | 0.857 [0.853, 0.862] | 0.827 [0.817, 0.836] | 0.850 [0.841, 0.858] |
| **swap: Null B** | 0.694 [0.669, 0.718] | 0.677 [0.668, 0.686] | 0.602 [0.557, 0.646] | 0.628 [0.590, 0.666] |
| swap: B2 | 0.733 [0.721, 0.746] | 0.779 [0.765, 0.792] | 0.692 [0.657, 0.727] | 0.712 [0.663, 0.761] |
| 7th/8th: Null A | 0.821 [0.810, 0.832] | 0.886 [0.885, 0.888] | 0.850 [0.839, 0.860] | 0.916 [0.906, 0.925] |
| **7th/8th: Null B** | 0.696 [0.624, 0.769] | 0.726 [0.710, 0.742] | 0.690 [0.673, 0.708] | 0.693 [0.685, 0.701] |
| 7th/8th: B2 | 0.767 [0.727, 0.807] | 0.814 [0.799, 0.828] | 0.726 [0.708, 0.743] | 0.794 [0.774, 0.814] |

## F2 power (planted edges, Tabula Sapiens Mammary; AMENDMENT_1.md)

Pairs from LIANA+'s consensus resource. Power in the top 1% / AUROC planted vs one-sided / precision at k (k = true rows), mean of the replicates run (planted rows present: f = 1.5: 30.0 of 30, f = 2.0: 30.0 of 30).

| Variant | f = 1.5 | f = 2.0 |
|---|---|---|
| LIANA+ as shipped (5 methods) | 0.96 / 0.928 / 0.35 | 0.94 / 0.930 / 0.46 |
| **+ CellChat (baseline)** | 1.00 / 0.947 / 0.18 | 1.00 / 0.940 / 0.25 |
| **swap: Null B** | 1.00 / 0.996 / 0.70 | 1.00 / 0.998 / 0.78 |
| swap: B2 | 0.97 / 0.969 / 0.53 | 0.99 / 0.981 / 0.64 |
| 7th/8th: Null A | 0.99 / 0.939 / 0.11 | 1.00 / 0.925 / 0.08 |
| **7th/8th: Null B** | 1.00 / 0.993 / 0.46 | 1.00 / 0.994 / 0.54 |
| 7th/8th: B2 | 0.90 / 0.944 / 0.35 | 0.91 / 0.953 / 0.43 |

## F3 biology

Odds ratio among the baseline's calls (s ≠ t), kept vs dropped by the variant (95% bootstrap CI; ↑ above 1, ↓ below 1; – where it keeps or drops none; k where it keeps fewer than 20).

| Variant | tnbc_wu2021 spatial | allen_brain_atlas spatial | moffitt2018_hypothalamus spatial | tnbc_wu2021 cytosig | tnbc_wu2021 progeny | allen_brain_atlas cytosig | allen_brain_atlas progeny | tabula_sapiens_mammary cytosig | tabula_sapiens_mammary progeny | tabula_sapiens_liver cytosig | tabula_sapiens_liver progeny | moffitt2018_hypothalamus cytosig | moffitt2018_hypothalamus progeny |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| LIANA+ as shipped (5 methods) | 0.88 [0.57, 1.53] | 1.58 [0.86, 3.47] | 1.31 [0.51, 2.66] | 1.11 [0.72, 2.28] | 2.64 [1.64, 7.44]↑ | 5.52 [2.73, 33.93]↑ | 1.23 [0.24, 6.00] | 8.34 [2.94, 38.16]↑ | 2.63 [0.68, 22.09] | 1.29 [0.90, 3.57] | 10.82 [3.93, 27.06]↑ | 2.17 [0.18, 6.18] | 0.11 [0.03, 1.21] (k=16) |
| **swap: Null B** | 0.50 [0.29, 1.44] | 1.35 [0.96, 1.86] | 1.15 [0.72, 2.13] | 0.51 [0.21, 1.40] | 6.28 [2.24, 29.91]↑ | 5.53 [2.80, 15.93]↑ | 0.63 [0.19, 5.41] | 1.51 [0.63, 16.75] | 0.77 [0.47, 1.52] | 0.79 [0.32, 2.64] | 0.47 [0.21, 0.95]↓ | – | 0.66 [0.02, 2.90] |
| swap: B2 | 1.33 [0.80, 2.22] | 1.48 [0.96, 2.39] | 1.82 [0.93, 4.33] | 0.44 [0.18, 2.44] | 4.43 [1.21, 29.01]↑ | 2.94 [1.50, 7.04]↑ | 0.47 [0.12, 3.44] | 1.16 [0.48, 10.54] | 0.71 [0.35, 2.91] | 0.61 [0.28, 7.92] | 0.42 [0.15, 2.66] | 0.16 [0.01, 0.51]↓ | 0.38 [0.02, 1.51] |
| 7th/8th: Null A | 0.95 [0.34, 6.02] | 0.51 [0.36, 0.85]↓ | 0.88 [0.44, 3.03] | 0.90 [0.27, 3.73] | 1.31 [0.26, 12.15] | 0.26 [0.15, 0.50]↓ | 0.97 [0.20, 3.69] | 0.27 [0.04, 1.13] | 1.88 [0.25, 5.46] | 1.13 [0.30, 7.65] | 0.39 [0.11, 3.47] | 0.02 [0.01, 0.23]↓ | – |
| **7th/8th: Null B** | 1.01 [0.48, 1.97] | 1.07 [0.81, 1.35] | 1.10 [0.72, 1.88] | 0.59 [0.23, 1.86] | 3.77 [1.15, 15.14]↑ | 1.83 [1.15, 3.34]↑ | 0.64 [0.19, 5.63] | 0.76 [0.35, 4.16] | 0.60 [0.35, 1.41] | 0.78 [0.27, 3.01] | 0.24 [0.09, 0.70]↓ | 0.02 [0.01, 0.23]↓ | 0.66 [0.02, 2.90] |
| 7th/8th: B2 | 1.15 [0.62, 1.95] | 1.24 [0.93, 1.64] | 1.36 [0.78, 2.81] | 0.95 [0.33, 5.87] | 3.96 [1.21, 25.44]↑ | 1.78 [1.13, 3.50]↑ | 0.98 [0.42, 5.82] | 0.71 [0.29, 6.32] | 0.39 [0.19, 1.92] | 1.14 [0.30, 11.30] | 0.74 [0.25, 4.39] | 0.02 [0.00, 0.19]↓ | 0.38 [0.02, 1.51] |

Pair-level AUPRC (pairs scored by their number of calls; difference from the baseline, paired bootstrap):

| Dataset | LIANA+ as shipped (5 methods) | **+ CellChat (baseline)** | **swap: Null B** | swap: B2 | 7th/8th: Null A | **7th/8th: Null B** | 7th/8th: B2 |
|---|---|---|---|---|---|---|---|
| tnbc_wu2021 | 0.015 (+0.001) | 0.014 | 0.014 (-0.001) | 0.015 (+0.000) | 0.015 (+0.001) | 0.015 (+0.000) | 0.014 (-0.000) |
| allen_brain_atlas | 0.056 (+0.005) | 0.051 | 0.059 (+0.008) | 0.059 (+0.007) | 0.043 (-0.008) | 0.050 (-0.001) | 0.049 (-0.002) |
| moffitt2018_hypothalamus | 0.249 (-0.006) | 0.255 | 0.261 (+0.006) | 0.293 (+0.039) | 0.236 (-0.018) | 0.243 (-0.012) | 0.281 (+0.026) |

## Decision (AMENDMENT_1.md's rule)

| Variant | Calibration | Power | Biology | Verdict |
|---|---|---|---|---|
| **swap: Null B** | yes | yes | no | valid, not biologically distinguishable |
| **7th/8th: Null B** | yes | yes | no | valid, not biologically distinguishable |
