# AMENDMENT_4 results: T2 on the rest of AMENDMENT_2

Plan: [`AMENDMENT_4.md`](../../../liana_final/AMENDMENT_4.md). T2 as in AMENDMENT_3; every run also produced liana6 and the Null B variants, so the columns are paired. Intervals: 95% t-intervals over replicates. Power cells: power in the top 1% / AUROC planted vs one-sided.

## Lung (A2.1): held-out confirmation

| Variant | 1 calibration | 2 power | Outcome |
|---|---|---|---|
| swap_null_b | holds | holds | confirmed |
| extra_null_b | holds | holds | confirmed |
| swap_t2 | fails | holds | not confirmed (calibration) |
| extra_t2 | fails | holds | not confirmed (calibration) |

| Scenario | reps | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|---|
| receptors shuffled | 3 | 0.880 [0.877, 0.883] | 0.677 [0.666, 0.687] | 0.753 [0.741, 0.764] | 0.873 [0.870, 0.876] | 0.885 [0.882, 0.888] |
| ligands shuffled | 3 | 0.833 [0.829, 0.837] | 0.679 [0.668, 0.691] | 0.733 [0.725, 0.741] | 0.892 [0.891, 0.893] | 0.894 [0.893, 0.896] |

| Effect | reps | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|---|
| f = 1.5 | 5 | 0.99 / 0.945 | 0.99 / 0.987 | 0.99 / 0.986 | 0.99 / 0.961 | 0.99 / 0.957 |
| f = 2.0 | 5 | 1.00 / 0.945 | 1.00 / 0.997 | 1.00 / 0.994 | 1.00 / 0.982 | 1.00 / 0.960 |

## Donors (A2.3 a): labels shuffled within donor, share of rows with p ≤ 0.05 (nominal 0.05)

| Data set | Tool | reps | Null A | Null B | T2 | T2 untestable |
|---|---|---|---|---|---|---|
| Lung | CellPhoneDB | 3 | 0.236 [0.222, 0.250] | 0.060 [0.051, 0.070] | 0.062 [0.059, 0.066] | 0.0% |
| Lung | CellChat | 3 | 0.297 [0.273, 0.321] | 0.075 [0.060, 0.091] | 0.125 [0.121, 0.130] | 0.0% |
| TNBC | CellPhoneDB | 3 | 0.369 [0.364, 0.374] | 0.128 [0.120, 0.136] | 0.057 [0.055, 0.060] | 0.0% |
| TNBC | CellChat | 3 | 0.506 [0.487, 0.524] | 0.184 [0.173, 0.196] | 0.228 [0.216, 0.240] | 0.0% |
| Mammary | CellPhoneDB | 3 | 0.178 [0.141, 0.216] | 0.037 [0.028, 0.046] | 0.055 [0.047, 0.064] | 0.0% |
| Mammary | CellChat | 3 | 0.193 [0.154, 0.232] | 0.040 [0.034, 0.046] | 0.134 [0.118, 0.150] | 0.0% |
| Liver | CellPhoneDB | 3 | 0.304 [0.279, 0.330] | 0.107 [0.097, 0.117] | 0.058 [0.052, 0.063] | 0.0% |
| Liver | CellChat | 3 | 0.377 [0.332, 0.422] | 0.146 [0.124, 0.169] | 0.123 [0.119, 0.127] | 0.0% |

## Donors (A2.3 b): real labels, within-donor vs global permutations

| Data set | swap_null_b: top-1% / call Jaccard | swap_t2: top-1% / call Jaccard | extra_t2: top-1% / call Jaccard |
|---|---|---|---|
| TNBC | 0.746 / 0.913 | 1.000 / 1.000 | 1.000 / 0.979 |
| Mammary | 0.799 / 0.913 | 1.000 / 1.000 | 0.892 / 0.961 |
| Liver | 0.597 / 0.868 | 1.000 / 1.000 | 0.853 / 0.863 |

## Sensitivity (A2.5)

Calibration: false share in the top 1% (below liana6?):

| Run | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|
| tabula_sapiens_mammary@expr_prop=0.05|receptor|0 | 0.849 | 0.744 ✓ | 0.726 ✓ | 0.840 ✓ | 0.852 ✗ |
| tabula_sapiens_mammary@expr_prop=0.05|ligand|0 | 0.874 | 0.768 ✓ | 0.812 ✓ | 0.917 ✗ | 0.917 ✗ |
| tabula_sapiens_mammary@expr_prop=0.2|receptor|0 | 0.734 | 0.566 ✓ | 0.535 ✓ | 0.814 ✗ | 0.783 ✗ |
| tabula_sapiens_mammary@expr_prop=0.2|ligand|0 | 0.752 | 0.509 ✓ | 0.555 ✓ | 0.852 ✗ | 0.915 ✗ |
| tabula_sapiens_mammary@cells=0.5|receptor|0 | 0.815 | 0.688 ✓ | 0.700 ✓ | 0.798 ✓ | 0.843 ✗ |
| tabula_sapiens_mammary@cells=0.25|receptor|0 | 0.771 | 0.655 ✓ | 0.713 ✓ | 0.778 ✗ | 0.808 ✗ |
| tabula_sapiens_mammary@labels=broad_cell_class|receptor|0 | 0.749 | 0.632 ✓ | 0.642 ✓ | 0.790 ✗ | 0.778 ✗ |
| tabula_sapiens_mammary@labels=broad_cell_class|ligand|0 | 0.853 | 0.719 ✓ | 0.748 ✓ | 0.892 ✗ | 0.901 ✗ |
| tabula_sapiens_mammary@labels=free_annotation|receptor|0 | 0.802 | 0.694 ✓ | 0.736 ✓ | 0.831 ✗ | 0.819 ✗ |
| tabula_sapiens_mammary@labels=free_annotation|ligand|0 | 0.839 | 0.684 ✓ | 0.722 ✓ | 0.898 ✗ | 0.904 ✗ |

Permutations, calibration (receptors, rep 0):

| Run | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|
| moffitt2018_hypothalamus@n_perms=1000|receptor|0 | 0.758 | 0.192 ✓ | 0.208 ✓ | 0.508 ✓ | 0.458 ✓ |
| moffitt2018_hypothalamus@n_perms=5000|receptor|0 | 0.758 | 0.183 ✓ | 0.225 ✓ | 0.508 ✓ | 0.458 ✓ |
| tabula_sapiens_mammary@n_perms=1000|receptor|0 | 0.804 | 0.694 ✓ | 0.673 ✓ | 0.817 ✗ | 0.817 ✗ |
| tabula_sapiens_mammary@n_perms=5000|receptor|0 | 0.806 | 0.695 ✓ | 0.653 ✓ | 0.817 ✗ | 0.838 ✗ |

Power, f = 1.5:

| Run | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|
| tabula_sapiens_mammary@expr_prop=0.05|f2|f1.5|0 | 1.00 / 0.961 | 1.00 / 0.997 | 1.00 / 0.995 | 1.00 / 0.987 | 0.97 / 0.970 |
| tabula_sapiens_mammary@expr_prop=0.05|f2|f1.5|1 | 1.00 / 0.956 | 1.00 / 0.997 | 1.00 / 0.992 | 1.00 / 0.977 | 1.00 / 0.961 |
| tabula_sapiens_mammary@expr_prop=0.2|f2|f1.5|0 | 1.00 / 0.933 | 1.00 / 0.995 | 1.00 / 0.994 | 0.90 / 0.956 | 0.77 / 0.897 |
| tabula_sapiens_mammary@expr_prop=0.2|f2|f1.5|1 | 0.90 / 0.914 | 1.00 / 0.994 | 1.00 / 0.986 | 0.87 / 0.943 | 0.70 / 0.885 |
| tabula_sapiens_mammary@cells=0.5|f2|f1.5|0 | 0.97 / 0.949 | 1.00 / 0.993 | 1.00 / 0.988 | 0.97 / 0.962 | 0.90 / 0.944 |
| tabula_sapiens_mammary@cells=0.25|f2|f1.5|0 | 0.90 / 0.920 | 0.93 / 0.977 | 0.93 / 0.974 | 0.87 / 0.927 | 0.83 / 0.922 |

Seeds, real labels, against seed 1337 (Spearman of specificity_rank / top-1% Jaccard):

| Run | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|
| seed1 | 0.9979 / 0.964 | 0.9965 / 0.935 | 0.9947 / 0.961 | 1.0000 / 1.000 | 0.9980 / 0.988 |
| seed2 | 0.9977 / 0.979 | 0.9960 / 0.973 | 0.9943 / 0.964 | 1.0000 / 1.000 | 0.9981 / 0.984 |

5,000 vs 1,000 permutations, real labels (Spearman of specificity_rank / top-1% Jaccard):

| Run | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|
| real | 0.9991 / 1.000 | 0.9969 / 1.000 | 0.9975 / 0.889 | 1.0000 / 1.000 | 0.9992 / 1.000 |
| real@tabula_sapiens_mammary | 0.9971 / 0.755 | 0.9963 / 0.734 | 0.9944 / 0.856 | 1.0000 / 1.000 | 0.9978 / 0.745 |

## Mouse (A2.7): Allen brain, criterion 1

| Variant | criterion 1 |
|---|---|
| swap_null_b | holds |
| extra_null_b | holds |
| swap_t2 | fails |
| extra_t2 | holds |

| Scenario | reps | liana6 | swap_null_b | extra_null_b | swap_t2 | extra_t2 |
|---|---|---|---|---|---|---|
| receptors shuffled | 3 | 0.836 [0.831, 0.840] | 0.726 [0.708, 0.744] | 0.690 [0.656, 0.723] | 0.874 [0.870, 0.879] | 0.831 [0.823, 0.839] |
| ligands shuffled | 3 | 0.925 [0.925, 0.926] | 0.825 [0.821, 0.829] | 0.808 [0.776, 0.839] | 0.939 [0.938, 0.939] | 0.902 [0.899, 0.905] |

## What swap_t2 drops (among liana6's calls, s ≠ t)

Breadth: share of cell types where the gene is detected in ≥ 10% of cells; a row's breadth is the larger of its ligand's and receptor's. p: one-sided Mann–Whitney, dropped broader than kept.

| Data set | kept | dropped | breadth dropped / kept | p | dropped by both / T2 only / Null B only | Jaccard of dropped sets | promoted |
|---|---|---|---|---|---|---|---|
| TNBC | 20481 | 1801 | 0.483 / 0.690 | 1 | 652 / 1149 / 3274 | 0.128 | 6034 |
| Allen brain | 36313 | 4586 | 0.826 / 0.913 | 1 | 2048 / 2538 / 7395 | 0.171 | 8644 |
| Mammary | 12687 | 1959 | 0.688 / 0.875 | 1 | 869 / 1090 / 2599 | 0.191 | 2159 |
| Liver | 28583 | 2790 | 0.917 / 0.958 | 1 | 785 / 2005 / 3609 | 0.123 | 7215 |
| hypothalamus | 1065 | 164 | 1.000 / 1.000 | 0.24 | 68 / 96 / 78 | 0.281 | 235 |

Tables: `cases_<data set>_dropped.csv` and `cases_<data set>_promoted.csv` in this folder.

## Timing

| Data set | T2 (s) | this pipeline, A2.10 (s) | rows (CellPhoneDB / CellChat) |
|---|---|---|---|
| Mammary | 2.02 | 2208.5 | 86,001 / 86,001 |
| hypothalamus | 0.75 | 1562.3 | 10,197 / 10,197 |

## Checks: liana6 and the Null B variants reproduce AMENDMENT_2

| File | Result |
|---|---|
| confirm_calibration | 6 of 6 runs identical |
| confirm_power | 10 of 10 runs identical |
| mouse_calibration | 6 of 6 runs identical |
| sensitivity_calibration | 10 of 10 runs identical |
| sensitivity_power | 6 of 6 runs identical |
| sensitivity_permutations_calibration | 4 of 4 runs identical |
| donor_null (Null A, Null B, global) | 12 of 12 runs identical |

