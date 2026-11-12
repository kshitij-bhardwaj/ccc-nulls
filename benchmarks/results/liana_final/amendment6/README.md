# AMENDMENT_6 results: Connectome's p-value replaced too

Plan: [`AMENDMENT_6.md`](../../../liana_final/AMENDMENT_6.md). `swap_<n>`: CellPhoneDB's and CellChat's p replaced by null n's; `swap3_<n>`: Connectome's specificity too. Same runs, so the columns are paired. Intervals: 95% t-intervals over replicates.

## Calibration: false share of the top 1% (lower is better)

| Setting | reps | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|---|
| Mammary, receptors | 3 | 0.799 [0.782, 0.815] | 0.694 [0.669, 0.718] | 0.733 [0.721, 0.746] | 0.818 [0.794, 0.843] | 0.813 [0.805, 0.822] | 0.499 [0.459, 0.539] | 0.694 [0.678, 0.710] | 0.832 [0.818, 0.847] |
| Mammary, ligands | 3 | 0.857 [0.853, 0.862] | 0.677 [0.668, 0.686] | 0.779 [0.765, 0.792] | 0.894 [0.884, 0.905] | 0.895 [0.889, 0.900] | 0.600 [0.579, 0.620] | 0.784 [0.774, 0.794] | 0.929 [0.924, 0.934] |
| Liver, receptors | 3 | 0.827 [0.817, 0.836] | 0.602 [0.557, 0.646] | 0.692 [0.657, 0.727] | 0.796 [0.771, 0.820] | 0.840 [0.833, 0.847] | 0.322 [0.300, 0.344] | 0.626 [0.568, 0.684] | 0.820 [0.788, 0.853] |
| Liver, ligands | 3 | 0.850 [0.841, 0.858] | 0.628 [0.590, 0.666] | 0.712 [0.663, 0.761] | 0.864 [0.860, 0.867] | 0.918 [0.912, 0.924] | 0.467 [0.439, 0.495] | 0.764 [0.730, 0.798] | 0.903 [0.897, 0.910] |
| Lung (held out), receptors | 3 | 0.880 [0.877, 0.883] | 0.677 [0.666, 0.687] | 0.758 [0.747, 0.769] | 0.873 [0.870, 0.876] | 0.902 [0.898, 0.906] | 0.434 [0.431, 0.438] | 0.701 [0.677, 0.725] | 0.875 [0.872, 0.878] |
| Lung (held out), ligands | 3 | 0.833 [0.829, 0.837] | 0.679 [0.668, 0.691] | 0.738 [0.733, 0.744] | 0.892 [0.891, 0.893] | 0.921 [0.918, 0.924] | 0.494 [0.485, 0.503] | 0.739 [0.730, 0.748] | 0.911 [0.910, 0.912] |
| Allen brain, receptors | 3 | 0.836 [0.831, 0.840] | 0.726 [0.708, 0.744] | 0.735 [0.725, 0.745] | 0.874 [0.870, 0.879] | 0.844 [0.839, 0.848] | 0.467 [0.459, 0.475] | 0.562 [0.547, 0.577] | 0.897 [0.893, 0.901] |
| Allen brain, ligands | 3 | 0.925 [0.925, 0.926] | 0.825 [0.821, 0.829] | 0.853 [0.849, 0.857] | 0.939 [0.938, 0.939] | 0.899 [0.896, 0.902] | 0.532 [0.523, 0.540] | 0.744 [0.716, 0.771] | 0.952 [0.949, 0.956] |

## Held-out Lung: AMENDMENT_1's rule

| Variant | 1 calibration | 2 power | Outcome |
|---|---|---|---|
| swap_null_b | holds | holds | confirmed |
| swap_b2 | holds | holds | confirmed |
| swap_t2 | fails | holds | not confirmed (calibration) |
| swap3_null_a | fails | fails | not confirmed (calibration, power) |
| swap3_null_b | holds | holds | confirmed |
| swap3_b2 | holds | fails | not confirmed (power) |
| swap3_t2 | fails | holds | not confirmed (calibration) |

## Power, F2 design: power in the top 1% / AUROC, single pairs

| Effect | reps | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|---|
| f = 1.5 | 5 | 1.00 / 0.947 | 1.00 / 0.996 | 0.97 / 0.969 | 0.99 / 0.982 | 0.99 / 0.933 | 1.00 / 0.999 | 0.75 / 0.923 | 0.94 / 0.956 |
| f = 2.0 | 5 | 1.00 / 0.940 | 1.00 / 0.998 | 0.99 / 0.981 | 0.99 / 0.985 | 1.00 / 0.915 | 1.00 / 1.000 | 0.80 / 0.946 | 0.95 / 0.967 |

## Power, small effects (A2.2): power in the top 1% / AUROC, single pairs

| Effect | reps | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|---|
| f = 1.1 | 5 | 0.55 / 0.868 | 0.65 / 0.914 | 0.54 / 0.874 | 0.46 / 0.872 | 0.50 / 0.859 | 0.71 / 0.926 | 0.50 / 0.821 | 0.43 / 0.860 |
| f = 1.25 | 5 | 0.73 / 0.915 | 0.99 / 0.985 | 0.83 / 0.944 | 0.89 / 0.952 | 0.61 / 0.894 | 0.99 / 0.991 | 0.72 / 0.898 | 0.85 / 0.937 |
| f = 1.5 | 5 | 0.98 / 0.951 | 1.00 / 0.996 | 0.93 / 0.972 | 0.97 / 0.978 | 0.97 / 0.937 | 1.00 / 0.999 | 0.79 / 0.934 | 0.93 / 0.965 |

## Power, held-out Lung: power in the top 1% / AUROC, single pairs

| Effect | reps | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|---|
| f = 1.5 | 5 | 0.99 / 0.945 | 0.99 / 0.987 | 0.92 / 0.951 | 0.99 / 0.961 | 0.99 / 0.942 | 0.99 / 0.991 | 0.70 / 0.905 | 0.97 / 0.951 |
| f = 2.0 | 5 | 1.00 / 0.945 | 1.00 / 0.997 | 0.99 / 0.972 | 1.00 / 0.982 | 1.00 / 0.939 | 1.00 / 0.998 | 0.82 / 0.938 | 0.98 / 0.970 |

## Connectome's own p on false rows (Mammary and Liver; share ≤ 0.05 / BH discoveries, false share)

| Null | runs | p ≤ 0.05 | BH discoveries | false share of discoveries |
|---|---|---|---|---|
| null_a | 12 | 0.393 [0.372, 0.415] | 72,864 | 0.940 |
| null_b | 12 | 0.026 [0.023, 0.029] | 0 | – |
| b2 | 12 | 0.001 [0.000, 0.002] | 2 | 0.000 |
| t2 | 12 | 0.073 [0.070, 0.075] | 1,509 | 0.908 |

## Real labels: calls (specificity_rank ≤ 0.05) and Jaccard with liana6

| Data set | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|
| tnbc_wu2021 | 23,466 (1.00) | 24,005 (0.69) | 34,458 (0.59) | 27,935 (0.72) | 21,885 (0.88) | 25,068 (0.61) | 37,853 (0.51) | 28,154 (0.64) |
| allen_brain_atlas | 42,879 (1.00) | 39,249 (0.67) | 50,911 (0.62) | 47,145 (0.74) | 40,580 (0.91) | 39,720 (0.57) | 56,037 (0.51) | 46,101 (0.67) |
| tabula_sapiens_mammary | 15,946 (1.00) | 13,208 (0.72) | 19,373 (0.64) | 16,151 (0.76) | 16,883 (0.89) | 14,284 (0.62) | 21,254 (0.54) | 16,605 (0.66) |
| tabula_sapiens_liver | 33,349 (1.00) | 33,740 (0.75) | 45,527 (0.65) | 37,991 (0.74) | 31,007 (0.90) | 33,799 (0.69) | 48,921 (0.56) | 37,669 (0.68) |
| moffitt2018_hypothalamus | 1,456 (1.00) | 1,292 (0.89) | 1,937 (0.65) | 1,545 (0.74) | 1,492 (0.80) | 1,287 (0.69) | 2,337 (0.48) | 1,670 (0.52) |

## Biology (A2.8): spatial odds ratio, kept vs dropped [95% CI], q across all

| Variant | TNBC | Allen brain | hypothalamus |
|---|---|---|---|
| swap_null_b | 0.50 [0.29, 1.44] q 0.47 | 1.35 [0.96, 1.86] q 0.34 | 1.15 [0.72, 2.13] q 0.83 |
| swap_b2 | 1.33 [0.80, 2.22] q 0.54 | 1.48 [0.96, 2.39] q 0.32 | 1.82 [0.93, 4.33] q 0.34 |
| swap_t2 | 0.88 [0.64, 1.27] q 0.77 | 1.15 [0.86, 1.61] q 0.74 | 1.03 [0.68, 1.53] q 0.90 |
| swap3_null_a | 1.09 [0.42, 2.81] q 0.90 | 0.55 [0.37, 0.95] q 0.21 | 0.61 [0.33, 1.88] q 0.59 |
| swap3_null_b | 0.54 [0.31, 1.75] q 0.54 | 1.41 [1.06, 1.80] q 0.11 | 0.89 [0.55, 1.71] q 0.86 |
| swap3_b2 | 1.33 [0.76, 2.09] q 0.60 | 1.22 [0.95, 1.54] q 0.40 | 1.13 [0.66, 2.78] q 0.87 |
| swap3_t2 | 1.15 [0.74, 1.80] q 0.77 | 0.99 [0.84, 1.21] q 0.90 | 0.72 [0.48, 1.23] q 0.53 |

Tests with q (all) < 0.05: swap_null_b | tnbc_wu2021 progeny, swap3_null_b | tnbc_wu2021 progeny, swap3_b2 | tnbc_wu2021 progeny, swap_null_b | allen_brain_atlas cytosig, swap_b2 | allen_brain_atlas cytosig, swap3_null_a | allen_brain_atlas cytosig, swap3_null_b | allen_brain_atlas cytosig, swap_t2 | tabula_sapiens_mammary cytosig, swap3_null_b | tabula_sapiens_liver progeny, swap_b2 | moffitt2018_hypothalamus cytosig, swap3_null_a | moffitt2018_hypothalamus cytosig, swap3_null_b | moffitt2018_hypothalamus cytosig, swap3_b2 | moffitt2018_hypothalamus cytosig.

Check: 34 of 34 tests that AMENDMENT_5 also ran are identical (odds ratio and p). The calls are joined to the downstream rows in human symbols, as `evaluate.py` built them; AMENDMENT_5 joined T2's calls in mouse symbols, so its `swap_t2` and `extra_t2` downstream tests on Allen brain and the hypothalamus are empty. Here `swap_t2` has them.

## Checks: liana6 and the earlier variants reproduce AMENDMENTS 3 and 4

| File | runs identical |
|---|---|
| calibration | 12 of 12 |
| confirm_calibration | 6 of 6 |
| mouse_calibration | 6 of 6 |
| power | 10 of 10 |
| hard_power | 15 of 15 |
| confirm_power | 10 of 10 |

