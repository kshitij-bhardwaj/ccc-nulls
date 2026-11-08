# AMENDMENT_3 results: T2 inside LIANA+'s consensus

T2 (PR #11): each row's score against the same interaction's scores at its other cell-type pairs (leave-one-out t). Plan: [`AMENDMENT_3.md`](../../../liana_final/AMENDMENT_3.md). Every number below comes from runs that also produced liana6 and the Null B variants, so the columns are paired. Intervals: 95% t-intervals over replicates. Where each run was computed: the `machine` field of each run in the JSON files.

## Decision (AMENDMENT_1's rule)

| Variant | 1 calibration | 2 power | 3 biology | Verdict |
|---|---|---|---|---|
| swap_t2 | False | True | False | not supported |
| extra_t2 | False | True | False | not supported |

For comparison, both Null B variants: *valid, not biologically distinguishable* (PLAN.md's results).

## Calibration: false share in the top 1% (lower is better)

| Scenario | reps | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|---|
| Mammary, receptors shuffled | 3 | 0.799 [0.782, 0.815] | 0.694 [0.669, 0.718] | 0.696 [0.624, 0.769] | 0.818 [0.794, 0.843] | 0.821 [0.811, 0.830] |
| Mammary, ligands shuffled | 3 | 0.857 [0.853, 0.862] | 0.677 [0.668, 0.686] | 0.726 [0.710, 0.742] | 0.894 [0.884, 0.905] | 0.922 [0.913, 0.930] |
| Liver, receptors shuffled | 3 | 0.827 [0.817, 0.836] | 0.602 [0.557, 0.646] | 0.690 [0.673, 0.708] | 0.796 [0.771, 0.820] | 0.829 [0.817, 0.840] |
| Liver, ligands shuffled | 3 | 0.850 [0.841, 0.858] | 0.628 [0.590, 0.666] | 0.693 [0.685, 0.701] | 0.864 [0.860, 0.867] | 0.895 [0.894, 0.897] |

## p-values of the false rows (share with p ≤ 0.05; a valid test gives ≤ 0.05)

| Tool | runs | p ≤ 0.05, null_a | p ≤ 0.05, null_b | p ≤ 0.05, b2 | p ≤ 0.05, t2 | BH 5% calls, null_a | BH 5% calls, null_b | BH 5% calls, b2 | BH 5% calls, t2 |
|---|---|---|---|---|---|---|---|---|---|
| CellPhoneDB | 12 | 0.383 [0.364, 0.402] | 0.026 [0.023, 0.029] | 0.001 [0.000, 0.002] | 0.073 [0.070, 0.077] | 70,597 | 0 | 2 | 423 |
| CellChat | 12 | 0.188 [0.167, 0.208] | 0.013 [0.011, 0.015] | 0.001 [0.000, 0.002] | 0.052 [0.050, 0.054] | 34,106 | 0 | 2 | 2,466 |

## Power (F2 design): power in the top 1% / AUROC planted vs one-sided

| Effect | pairs | reps | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|---|---|
| f = 1.5 | single | 5 | 1.00 / 0.947 | 1.00 / 0.996 | 1.00 / 0.993 | 0.99 / 0.982 | 0.95 / 0.955 |
| f = 2.0 | single | 5 | 1.00 / 0.940 | 1.00 / 0.998 | 1.00 / 0.994 | 0.99 / 0.985 | 0.95 / 0.955 |

## Hard power (A2.2 design): power in the top 1% / AUROC planted vs one-sided

| Effect | pairs | reps | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|---|---|
| f = 1.1 | single | 5 | 0.55 / 0.868 | 0.65 / 0.914 | 0.59 / 0.898 | 0.46 / 0.872 | 0.40 / 0.858 |
| f = 1.1 | complex | 5 | 0.56 / 0.864 | 0.62 / 0.912 | 0.60 / 0.898 | 0.34 / 0.850 | 0.30 / 0.835 |
| f = 1.1 | unequal | 5 | 0.30 / 0.722 | 0.36 / 0.779 | 0.36 / 0.758 | 0.24 / 0.731 | 0.26 / 0.706 |
| f = 1.25 | single | 5 | 0.73 / 0.915 | 0.99 / 0.985 | 0.89 / 0.960 | 0.89 / 0.952 | 0.67 / 0.897 |
| f = 1.25 | complex | 5 | 0.70 / 0.879 | 0.98 / 0.974 | 0.88 / 0.937 | 0.78 / 0.895 | 0.54 / 0.824 |
| f = 1.25 | unequal | 5 | 0.38 / 0.790 | 0.86 / 0.932 | 0.56 / 0.881 | 0.54 / 0.820 | 0.36 / 0.740 |
| f = 1.5 | single | 5 | 0.98 / 0.951 | 1.00 / 0.996 | 0.99 / 0.992 | 0.97 / 0.978 | 0.93 / 0.961 |
| f = 1.5 | complex | 5 | 0.96 / 0.946 | 0.98 / 0.996 | 1.00 / 0.992 | 0.80 / 0.934 | 0.78 / 0.904 |
| f = 1.5 | unequal | 5 | 0.62 / 0.838 | 0.98 / 0.976 | 0.80 / 0.946 | 0.80 / 0.905 | 0.60 / 0.827 |

## T2 − Null B on the same runs (reported, not judged)

Calibration: lower is better, so a negative difference favours T2. Power and AUROC: a positive difference favours T2.

| Difference | Measure | Mean [95% CI] |
|---|---|---|
| swap_t2 − swap_null_b | false share top 1%, Mammary receptor | 0.125 [0.122, 0.127] |
| swap_t2 − swap_null_b | false share top 1%, Mammary ligand | 0.217 [0.215, 0.219] |
| swap_t2 − swap_null_b | false share top 1%, Liver receptor | 0.194 [0.168, 0.220] |
| swap_t2 − swap_null_b | false share top 1%, Liver ligand | 0.236 [0.201, 0.271] |
| swap_t2 − swap_null_b | power_top, F2 single, f = 1.5 | -0.007 [-0.025, 0.012] |
| swap_t2 − swap_null_b | auroc, F2 single, f = 1.5 | -0.014 [-0.024, -0.004] |
| swap_t2 − swap_null_b | power_top, F2 single, f = 2.0 | -0.013 [-0.036, 0.009] |
| swap_t2 − swap_null_b | auroc, F2 single, f = 2.0 | -0.013 [-0.029, 0.002] |
| swap_t2 − swap_null_b | power_top, A2.2 single, f = 1.1 | -0.193 [-0.273, -0.114] |
| swap_t2 − swap_null_b | auroc, A2.2 single, f = 1.1 | -0.042 [-0.053, -0.031] |
| swap_t2 − swap_null_b | power_top, A2.2 single, f = 1.25 | -0.107 [-0.210, -0.004] |
| swap_t2 − swap_null_b | auroc, A2.2 single, f = 1.25 | -0.033 [-0.052, -0.014] |
| swap_t2 − swap_null_b | power_top, A2.2 single, f = 1.5 | -0.027 [-0.081, 0.027] |
| swap_t2 − swap_null_b | auroc, A2.2 single, f = 1.5 | -0.018 [-0.044, 0.009] |
| extra_t2 − extra_null_b | false share top 1%, Mammary receptor | 0.124 [0.061, 0.188] |
| extra_t2 − extra_null_b | false share top 1%, Mammary ligand | 0.196 [0.171, 0.220] |
| extra_t2 − extra_null_b | false share top 1%, Liver receptor | 0.138 [0.132, 0.144] |
| extra_t2 − extra_null_b | false share top 1%, Liver ligand | 0.202 [0.196, 0.209] |
| extra_t2 − extra_null_b | power_top, F2 single, f = 1.5 | -0.053 [-0.101, -0.006] |
| extra_t2 − extra_null_b | auroc, F2 single, f = 1.5 | -0.038 [-0.058, -0.019] |
| extra_t2 − extra_null_b | power_top, F2 single, f = 2.0 | -0.047 [-0.094, 0.001] |
| extra_t2 − extra_null_b | auroc, F2 single, f = 2.0 | -0.039 [-0.057, -0.020] |
| extra_t2 − extra_null_b | power_top, A2.2 single, f = 1.1 | -0.193 [-0.261, -0.125] |
| extra_t2 − extra_null_b | auroc, A2.2 single, f = 1.1 | -0.040 [-0.050, -0.030] |
| extra_t2 − extra_null_b | power_top, A2.2 single, f = 1.25 | -0.220 [-0.306, -0.134] |
| extra_t2 − extra_null_b | auroc, A2.2 single, f = 1.25 | -0.063 [-0.080, -0.045] |
| extra_t2 − extra_null_b | power_top, A2.2 single, f = 1.5 | -0.060 [-0.114, -0.006] |
| extra_t2 − extra_null_b | auroc, A2.2 single, f = 1.5 | -0.032 [-0.054, -0.010] |

## The five data sets: calls (specificity_rank ≤ 0.05) and Jaccard with liana6

| Data set | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 | T2 untestable rows (CellPhoneDB / CellChat) |
|---|---|---|---|---|---|---|
| TNBC | 23,466 (1.00) | 24,005 (0.69) | 23,122 (0.74) | 27,935 (0.72) | 27,946 (0.73) | 0.1% / 0.1% |
| Allen brain | 42,879 (1.00) | 39,249 (0.67) | 42,186 (0.68) | 47,145 (0.74) | 50,343 (0.73) | 0.2% / 0.2% |
| Mammary | 15,946 (1.00) | 13,208 (0.72) | 15,041 (0.73) | 16,151 (0.76) | 17,450 (0.75) | 0.4% / 0.4% |
| Liver | 33,349 (1.00) | 33,740 (0.75) | 34,632 (0.75) | 37,991 (0.74) | 38,757 (0.75) | 0.3% / 0.3% |
| hypothalamus | 1,456 (1.00) | 1,292 (0.89) | 1,480 (0.69) | 1,545 (0.74) | 1,600 (0.67) | 5.4% / 5.4% |

## Biology (F3): spatial odds ratio, liana6 calls kept vs dropped [95% CI]

| Data set | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|
| TNBC | 0.50 [0.29, 1.44] | 1.01 [0.48, 1.97] | 0.88 [0.64, 1.27] | 0.77 [0.43, 2.46] |
| Allen brain | 1.35 [0.96, 1.86] | 1.07 [0.81, 1.35] | 1.15 [0.86, 1.61] | 0.85 [0.69, 1.10] |
| hypothalamus | 1.15 [0.72, 2.13] | 1.10 [0.72, 1.88] | 1.03 [0.68, 1.53] | 0.90 [0.60, 1.78] |

Downstream response, odds ratio kept vs dropped:

| Data set | footprints | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|
| TNBC | cytosig | 0.515 | 0.588 | 0.941 | 0.531 |
| TNBC | progeny | 6.283 | 3.766 | 0.498 | 2.125 |
| Allen brain | cytosig | 5.531 | 1.835 | 1.790 | 1.189 |
| Allen brain | progeny | 0.630 | 0.636 | 0.498 | 0.352 |
| Mammary | cytosig | 1.511 | 0.761 | 3.549 | 1.050 |
| Mammary | progeny | 0.770 | 0.596 | 1.531 | 1.316 |
| Liver | cytosig | 0.786 | 0.784 | 1.259 | 2.502 |
| Liver | progeny | 0.472 | 0.241 | 0.787 | 0.265 |
| hypothalamus | cytosig | – | 0.019 | 1.000 | 0.235 |
| hypothalamus | progeny | 0.660 | 0.660 | – | 1.286 |

## Checks: the earlier variants reproduce

| Check | Result |
|---|---|
| calibration vs A2.6 | 12 of 12 runs identical |
| hard power vs A2.2 | 15 of 15 runs identical |
| power vs F2 (means) | largest difference 0.0000 |
| real TNBC vs run tables | {'liana6': True, 'swap_null_b': True, 'extra_null_b': True} |
| real Allen brain vs run tables | {'liana6': True, 'swap_null_b': True, 'extra_null_b': True} |
| real Mammary vs run tables | {'liana6': True, 'swap_null_b': True, 'extra_null_b': True} |
| real Liver vs run tables | {'liana6': True, 'swap_null_b': True, 'extra_null_b': True} |
| real hypothalamus vs run tables | {'liana6': True, 'swap_null_b': True, 'extra_null_b': True} |

