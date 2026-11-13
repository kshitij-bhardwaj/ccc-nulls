# AMENDMENT_7 results: Connectome through the remaining tests

Plan: [`AMENDMENT_7.md`](../../../liana_final/AMENDMENT_7.md). Reported, not decided: AMENDMENT_6's decision for `swap3_<n>` stands. Names as in the paper: JLP = Null A (LIANA+'s own test), BPT = Null B, MDT = B2, LOO-t = T2. Intervals: 95% t-intervals over replicates.

## Donors (A2.3 a): labels shuffled within donor; share of rows with p ≤ 0.05

s ≠ t, magnitude > 0. Global: labels permuted over all cells (LIANA+'s scheme); donor: within each donor. MDT uses no permutations, so it has the global column only.

| Data set | Tool | reps (global / donor) | JLP: global / donor | BPT: global / donor | LOO-t: global / donor | MDT |
|---|---|---|---|---|---|---|
| Lung | CellPhoneDB | 3 / 3 | 0.236 [0.222, 0.250] / 0.064 [0.051, 0.077] | 0.060 [0.051, 0.070] / 0.004 [0.002, 0.005] | 0.062 [0.059, 0.066] / 0.062 [0.059, 0.066] | 0.006 [0.005, 0.007] |
| Lung | CellChat | 3 / 3 | 0.297 [0.273, 0.321] / 0.073 [0.053, 0.093] | 0.075 [0.060, 0.091] / 0.004 [0.003, 0.005] | 0.125 [0.121, 0.130] / 0.125 [0.121, 0.130] | 0.010 [0.009, 0.012] |
| Lung | Connectome | 3 / 3 | 0.234 [0.222, 0.246] / 0.067 [0.052, 0.082] | 0.060 [0.051, 0.070] / 0.004 [0.002, 0.005] | 0.063 [0.061, 0.065] / 0.063 [0.061, 0.065] | 0.006 [0.005, 0.007] |
| TNBC | CellPhoneDB | 3 / 3 | 0.369 [0.364, 0.374] / 0.062 [0.047, 0.076] | 0.128 [0.120, 0.136] / 0.003 [0.002, 0.005] | 0.057 [0.055, 0.060] / 0.057 [0.055, 0.060] | 0.019 [0.016, 0.023] |
| TNBC | CellChat | 3 / 3 | 0.506 [0.487, 0.524] / 0.060 [0.041, 0.078] | 0.184 [0.173, 0.196] / 0.003 [0.002, 0.004] | 0.228 [0.216, 0.240] / 0.228 [0.216, 0.240] | 0.046 [0.045, 0.047] |
| TNBC | Connectome | 3 / 3 | 0.383 [0.369, 0.396] / 0.066 [0.053, 0.080] | 0.128 [0.120, 0.136] / 0.003 [0.002, 0.005] | 0.070 [0.069, 0.071] / 0.070 [0.069, 0.071] | 0.019 [0.016, 0.023] |
| Mammary | CellPhoneDB | 3 / 3 | 0.178 [0.141, 0.216] / 0.075 [0.056, 0.094] | 0.037 [0.028, 0.046] / 0.004 [0.003, 0.005] | 0.055 [0.047, 0.064] / 0.055 [0.047, 0.064] | 0.003 [0.002, 0.005] |
| Mammary | CellChat | 3 / 3 | 0.193 [0.154, 0.232] / 0.102 [0.073, 0.130] | 0.040 [0.034, 0.046] / 0.006 [0.004, 0.007] | 0.134 [0.118, 0.150] / 0.134 [0.118, 0.150] | 0.004 [0.001, 0.007] |
| Mammary | Connectome | 3 / 3 | 0.178 [0.138, 0.218] / 0.079 [0.060, 0.099] | 0.037 [0.028, 0.046] / 0.004 [0.003, 0.005] | 0.064 [0.058, 0.069] / 0.064 [0.058, 0.069] | 0.003 [0.002, 0.005] |
| Liver | CellPhoneDB | 3 / 3 | 0.304 [0.279, 0.330] / 0.060 [0.047, 0.073] | 0.107 [0.097, 0.117] / 0.003 [0.002, 0.003] | 0.058 [0.052, 0.063] / 0.058 [0.052, 0.063] | 0.004 [0.003, 0.004] |
| Liver | CellChat | 3 / 3 | 0.377 [0.332, 0.422] / 0.075 [0.051, 0.100] | 0.146 [0.124, 0.169] / 0.004 [0.002, 0.005] | 0.123 [0.119, 0.127] / 0.123 [0.119, 0.127] | 0.005 [0.003, 0.006] |
| Liver | Connectome | 3 / 3 | 0.286 [0.258, 0.313] / 0.064 [0.053, 0.075] | 0.107 [0.097, 0.117] / 0.003 [0.002, 0.003] | 0.062 [0.060, 0.065] / 0.062 [0.060, 0.065] | 0.004 [0.003, 0.004] |

## Sensitivity (A2.5), with the swap3 variants

✓: false share of the top 1% below liana6's in the same run; ✗: not below.

Calibration: false share in the top 1%:

| Run | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|
| tabula_sapiens_mammary@cells=0.5|receptor|0 | 0.815 | 0.688 ✓ | 0.730 ✓ | 0.798 ✓ | 0.833 ✗ | 0.510 ✓ | 0.700 ✓ | 0.810 ✓ |
| tabula_sapiens_mammary@cells=0.25|receptor|0 | 0.771 | 0.655 ✓ | 0.695 ✓ | 0.778 ✗ | 0.817 ✗ | 0.438 ✓ | 0.686 ✓ | 0.796 ✗ |
| tabula_sapiens_mammary@labels=broad_cell_class|receptor|0 | 0.749 | 0.632 ✓ | 0.656 ✓ | 0.790 ✗ | 0.794 ✗ | 0.311 ✓ | 0.602 ✓ | 0.813 ✗ |
| tabula_sapiens_mammary@labels=broad_cell_class|ligand|0 | 0.853 | 0.719 ✓ | 0.809 ✓ | 0.892 ✗ | 0.914 ✗ | 0.580 ✓ | 0.809 ✓ | 0.923 ✗ |
| tabula_sapiens_mammary@labels=free_annotation|receptor|0 | 0.802 | 0.694 ✓ | 0.731 ✓ | 0.831 ✗ | 0.828 ✗ | 0.523 ✓ | 0.710 ✓ | 0.843 ✗ |
| tabula_sapiens_mammary@labels=free_annotation|ligand|0 | 0.839 | 0.684 ✓ | 0.763 ✓ | 0.898 ✗ | 0.889 ✗ | 0.589 ✓ | 0.778 ✓ | 0.942 ✗ |
| tabula_sapiens_mammary@expr_prop=0.05|receptor|0 | 0.849 | 0.744 ✓ | 0.789 ✓ | 0.840 ✓ | 0.846 ✓ | 0.631 ✓ | 0.761 ✓ | 0.846 ✓ |
| tabula_sapiens_mammary@expr_prop=0.05|ligand|0 | 0.874 | 0.768 ✓ | 0.844 ✓ | 0.917 ✗ | 0.903 ✗ | 0.719 ✓ | 0.852 ✓ | 0.939 ✗ |
| tabula_sapiens_mammary@expr_prop=0.2|receptor|0 | 0.734 | 0.566 ✓ | 0.549 ✓ | 0.814 ✗ | 0.781 ✗ | 0.489 ✓ | 0.496 ✓ | 0.830 ✗ |
| tabula_sapiens_mammary@expr_prop=0.2|ligand|0 | 0.752 | 0.509 ✓ | 0.649 ✓ | 0.852 ✗ | 0.852 ✗ | 0.387 ✓ | 0.638 ✓ | 0.925 ✗ |

Permutations, calibration (receptors, rep 0):

| Run | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|
| moffitt2018_hypothalamus@n_perms=1000|receptor|0 | 0.758 | 0.192 ✓ | 0.442 ✓ | 0.508 ✓ | 0.625 ✓ | 0.117 ✓ | 0.250 ✓ | 0.533 ✓ |
| moffitt2018_hypothalamus@n_perms=5000|receptor|0 | 0.758 | 0.183 ✓ | 0.442 ✓ | 0.508 ✓ | 0.642 ✓ | 0.108 ✓ | 0.250 ✓ | 0.533 ✓ |
| tabula_sapiens_mammary@n_perms=1000|receptor|0 | 0.804 | 0.694 ✓ | 0.734 ✓ | 0.817 ✗ | 0.813 ✗ | 0.496 ✓ | 0.687 ✓ | 0.833 ✗ |
| tabula_sapiens_mammary@n_perms=5000|receptor|0 | 0.806 | 0.695 ✓ | 0.734 ✓ | 0.817 ✗ | 0.818 ✗ | 0.503 ✓ | 0.687 ✓ | 0.833 ✗ |

Calibration: runs where the variant's false share is below liana6's — swap_null_b 10 of 10, swap_b2 10 of 10, swap_t2 2 of 10, swap3_null_a 1 of 10, swap3_null_b 10 of 10, swap3_b2 10 of 10, swap3_t2 2 of 10.

Permutations, calibration: runs where the variant's false share is below liana6's — swap_null_b 4 of 4, swap_b2 4 of 4, swap_t2 2 of 4, swap3_null_a 2 of 4, swap3_null_b 4 of 4, swap3_b2 4 of 4, swap3_t2 2 of 4.

Power, f = 1.5 (power in the top 1% / AUROC, single pairs):

| Run | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|
| tabula_sapiens_mammary@cells=0.5|f2|f1.5|0 | 0.97 / 0.949 | 1.00 / 0.993 | 1.00 / 0.970 | 0.97 / 0.962 | 0.97 / 0.941 | 1.00 / 0.996 | 0.73 / 0.920 | 0.90 / 0.932 |
| tabula_sapiens_mammary@cells=0.25|f2|f1.5|0 | 0.90 / 0.920 | 0.93 / 0.977 | 0.93 / 0.946 | 0.87 / 0.927 | 0.93 / 0.922 | 0.97 / 0.989 | 0.70 / 0.908 | 0.83 / 0.921 |
| tabula_sapiens_mammary@expr_prop=0.05|f2|f1.5|0 | 1.00 / 0.961 | 1.00 / 0.997 | 1.00 / 0.963 | 1.00 / 0.987 | 1.00 / 0.951 | 1.00 / 0.999 | 0.63 / 0.906 | 0.97 / 0.973 |
| tabula_sapiens_mammary@expr_prop=0.05|f2|f1.5|1 | 1.00 / 0.956 | 1.00 / 0.997 | 0.97 / 0.976 | 1.00 / 0.977 | 1.00 / 0.936 | 1.00 / 0.999 | 0.80 / 0.928 | 0.93 / 0.963 |
| tabula_sapiens_mammary@expr_prop=0.2|f2|f1.5|0 | 1.00 / 0.933 | 1.00 / 0.995 | 0.83 / 0.944 | 0.90 / 0.956 | 0.93 / 0.913 | 1.00 / 0.999 | 0.63 / 0.879 | 0.77 / 0.889 |
| tabula_sapiens_mammary@expr_prop=0.2|f2|f1.5|1 | 0.90 / 0.914 | 1.00 / 0.994 | 0.93 / 0.965 | 0.87 / 0.943 | 0.87 / 0.874 | 1.00 / 0.999 | 0.67 / 0.900 | 0.60 / 0.862 |

Seeds, real labels (hypothalamus), against seed 1337 (Spearman of specificity_rank / top-1% Jaccard):

| Run | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|
| seed1 | 0.9979 / 0.964 | 0.9965 / 0.935 | not run | 1.0000 / 1.000 | 0.9968 / 0.970 | 0.9936 / 0.957 | not run | 1.0000 / 1.000 |
| seed2 | 0.9977 / 0.979 | 0.9960 / 0.973 | not run | 1.0000 / 1.000 | 0.9965 / 0.946 | 0.9931 / 0.977 | not run | 1.0000 / 1.000 |

5,000 vs 1,000 permutations, real labels (Spearman of specificity_rank / top-1% Jaccard):

| Run | liana6 | swap_null_b | swap_b2 | swap_t2 | swap3_null_a | swap3_null_b | swap3_b2 | swap3_t2 |
|---|---|---|---|---|---|---|---|---|
| real | 0.9991 / 1.000 | 0.9969 / 1.000 | not run | 1.0000 / 1.000 | 0.9982 / 0.889 | 0.9961 / 0.943 | not run | 1.0000 / 1.000 |
| real@tabula_sapiens_mammary | 0.9971 / 0.755 | 0.9963 / 0.734 | not run | 1.0000 / 1.000 | 0.9963 / 0.815 | 0.9953 / 0.792 | not run | 1.0000 / 1.000 |

"not run": A2.5's seed and permutation runs compute no MDT (`b2=False`, as in AMENDMENT_4), so `swap_b2` and `swap3_b2` have no second table to compare. The first seed run stopped on this (`KeyError`); `stability` was changed to record these variants as skipped (commit 42b952f), and the run was repeated. No other result changed: the sessions that had finished ran the same code for every part they computed.

## Case studies (A2.9): what swap3_<n> drops among liana6's calls (s ≠ t)

Breadth: share of cell types where the gene is detected in ≥ 10% of cells; a row's breadth is the larger of its ligand's and receptor's. p: one-sided Mann–Whitney, dropped broader than kept. Overlap with what swap_<n> drops in the same table. Computed locally from AMENDMENT_6's saved real tables, not in a Kaggle session (`python -m benchmarks.liana_final.amendment7 --part cases`).

| Data set | Variant | kept | dropped | breadth dropped / kept | p | dropped by both / swap3 only / swap only |
|---|---|---|---|---|---|---|
| TNBC | swap3_null_b | 17412 | 4870 | 0.897 / 0.586 | 0 | 3926 / 944 / 0 |
| TNBC | swap3_b2 | 19608 | 2674 | 0.897 / 0.586 | 2.9e-121 | 1961 / 713 / 0 |
| TNBC | swap3_t2 | 19118 | 3164 | 0.414 / 0.724 | 1 | 1781 / 1383 / 20 |
| Allen brain | swap3_null_b | 28480 | 12419 | 0.957 / 0.870 | 1.2e-146 | 9409 / 3010 / 34 |
| Allen brain | swap3_b2 | 31877 | 9022 | 0.913 / 0.870 | 2.9e-20 | 6895 / 2127 / 0 |
| Allen brain | swap3_t2 | 33993 | 6906 | 0.783 / 0.913 | 1 | 4571 / 2335 / 15 |
| Mammary | swap3_null_b | 10596 | 4050 | 0.875 / 0.812 | 2.3e-17 | 3458 / 592 / 10 |
| Mammary | swap3_b2 | 11959 | 2687 | 0.812 / 0.812 | 1 | 2010 / 677 / 17 |
| Mammary | swap3_t2 | 11786 | 2860 | 0.625 / 0.875 | 1 | 1948 / 912 / 11 |
| Liver | swap3_null_b | 25818 | 5555 | 0.958 / 0.958 | 7.5e-19 | 4389 / 1166 / 5 |
| Liver | swap3_b2 | 27840 | 3533 | 0.958 / 0.958 | 1 | 2129 / 1404 / 0 |
| Liver | swap3_t2 | 26979 | 4394 | 0.792 / 0.958 | 1 | 2784 / 1610 / 6 |
| hypothalamus | swap3_null_b | 938 | 291 | 1.000 / 0.875 | 0.68 | 146 / 145 / 0 |
| hypothalamus | swap3_b2 | 1044 | 185 | 0.875 / 1.000 | 0.96 | 105 / 80 / 0 |
| hypothalamus | swap3_t2 | 905 | 324 | 0.625 / 1.000 | 1 | 164 / 160 / 0 |

## Benchmark 1's tests at tool level, Connectome beside CellPhoneDB and CellChat

Each tool's own rows inside LIANA+, p from each null, BH at 5% over the rows with a p. JLP twice: LIANA+'s p (#/B, can be 0) and Benchmark 1's ((1 + #)/(B + 1)).

### Labels shuffled (every label permuted; any call is a false call)

Runs with any call / mean number of calls.

| Data set | Tool | reps | JLP (#/B) | JLP ((1+#)/(B+1)) | BPT | MDT | LOO-t |
|---|---|---|---|---|---|---|---|
| Mammary | CellPhoneDB | 3 | 3 / 164.3 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 758.7 |
| Mammary | CellChat | 3 | 3 / 64.3 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 3,289.3 |
| Mammary | Connectome | 3 | 3 / 179.7 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 2,349.0 |
| Liver | CellPhoneDB | 3 | 3 / 392.7 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 2,805.3 |
| Liver | CellChat | 3 | 3 / 222.3 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 4,449.0 |
| Liver | Connectome | 3 | 3 / 458.0 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 3,606.0 |

### Power (F2's design, 30 single pairs, f = 1.5, Mammary): power / FDP, mean over replicates

| B | Tool | reps | JLP (#/B) | JLP ((1+#)/(B+1)) | BPT | MDT | LOO-t |
|---|---|---|---|---|---|---|---|
| 1,000 | CellPhoneDB | 5 | 1.00 / 0.99 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.75 / 0.18 |
| 1,000 | CellChat | 5 | 0.99 / 0.97 | 0.00 / 0.00 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.86 |
| 1,000 | Connectome | 5 | 1.00 / 0.98 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.87 |
| 5,000 | CellPhoneDB | 5 | 1.00 / 0.98 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.75 / 0.18 |
| 5,000 | CellChat | 5 | 0.99 / 0.97 | 0.99 / 0.97 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.86 |
| 5,000 | Connectome | 5 | 1.00 / 0.98 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.87 |

### Real labels: BH calls

| Data set | Tool | rows | JLP (#/B) | JLP ((1+#)/(B+1)) | BPT | MDT | LOO-t |
|---|---|---|---|---|---|---|---|
| TNBC | CellPhoneDB | 154,213 | 106,730 | 106,646 | 61,810 | 357 | 11 |
| TNBC | CellChat | 154,213 | 40,856 | 40,817 | 25,448 | 357 | 3,296 |
| TNBC | Connectome | 154,213 | 113,459 | 113,375 | 61,810 | 357 | 1,886 |
| Allen brain | CellPhoneDB | 235,366 | 106,793 | 106,524 | 35,654 | 12 | 906 |
| Allen brain | CellChat | 235,366 | 69,719 | 69,570 | 30,557 | 12 | 4,916 |
| Allen brain | Connectome | 235,366 | 100,902 | 100,645 | 35,654 | 12 | 4,538 |
| Mammary | CellPhoneDB | 86,001 | 35,054 | 34,956 | 14,983 | 158 | 3 |
| Mammary | CellChat | 86,001 | 20,213 | 20,166 | 9,186 | 158 | 2,108 |
| Mammary | Connectome | 86,001 | 35,980 | 35,903 | 14,983 | 158 | 1,618 |
| Liver | CellPhoneDB | 196,381 | 88,904 | 88,704 | 39,270 | 9 | 296 |
| Liver | CellChat | 196,381 | 53,122 | 52,995 | 25,979 | 9 | 4,311 |
| Liver | Connectome | 196,381 | 90,109 | 89,897 | 39,270 | 9 | 3,255 |
| hypothalamus | CellPhoneDB | 10,197 | 6,565 | 6,562 | 4,232 | 10 | 5 |
| hypothalamus | CellChat | 10,197 | 2,645 | 2,643 | 1,838 | 10 | 146 |
| hypothalamus | Connectome | 10,197 | 6,925 | 6,923 | 4,232 | 10 | 104 |

### Spatial truth: among JLP's calls ((1 + #)/(B + 1)), odds ratio of the rows each null keeps vs drops [95% CI]; q: BH over all tests below

| Data set | Tool | BPT | MDT | LOO-t |
|---|---|---|---|---|
| TNBC | CellPhoneDB | 0.80 [0.60, 1.11] q 0.35 | 0.41 [0.09, 1.12] q 0.26 | 4.30 [1.93, 18.07] q 0.01 |
| TNBC | CellChat | 0.76 [0.59, 1.01] q 0.26 | 0.41 [0.09, 1.23] q 0.29 | 1.12 [0.33, 1.93] q 0.90 |
| TNBC | Connectome | 0.81 [0.66, 1.04] q 0.29 | 0.41 [0.09, 1.11] q 0.26 | 0.97 [0.23, 2.02] q 0.90 |
| Allen brain | CellPhoneDB | 1.39 [0.98, 1.90] q 0.26 | 1.92 [0.63, 7.97] q 0.50 | 2.36 [1.32, 3.72] q 0.08 |
| Allen brain | CellChat | 1.18 [0.96, 1.44] q 0.29 | 2.22 [0.72, 8.64] q 0.42 | 1.55 [0.90, 2.40] q 0.29 |
| Allen brain | Connectome | 1.35 [0.97, 1.80] q 0.26 | 1.88 [0.60, 7.86] q 0.52 | 2.05 [1.20, 2.98] q 0.10 |
| hypothalamus | CellPhoneDB | 0.89 [0.65, 1.21] q 0.58 | 3.09 [0.39, 19.74] q 0.50 | 2.05 [0.41, 20.29] q 0.53 |
| hypothalamus | CellChat | 1.04 [0.75, 1.48] q 0.90 | 3.31 [0.41, 20.25] q 0.49 | 1.03 [0.57, 1.77] q 0.90 |
| hypothalamus | Connectome | 0.84 [0.62, 1.15] q 0.46 | 3.00 [0.38, 19.18] q 0.50 | 1.27 [0.54, 2.42] q 0.61 |

Tests with q < 0.05: TNBC CellPhoneDB LOO-t.

## Checks: CellPhoneDB's and CellChat's numbers, and the earlier variants, reproduce earlier amendments

| Against | Result |
|---|---|
| AMENDMENT_2 (JLP, BPT; ≤ 0.01 and ≤ 0.05) | 24 of 24 runs identical |
| AMENDMENT_4 (LOO-t) | 12 of 12 runs identical |
| AMENDMENT_5 (MDT) | 12 of 12 runs identical |
| AMENDMENT_4 sensitivity_calibration (liana6, swap_null_b, swap_b2, swap_t2) | 10 of 10 runs identical |
| AMENDMENT_4 sensitivity_power (liana6, swap_null_b, swap_b2, swap_t2) | 6 of 6 runs identical |
| AMENDMENT_4 sensitivity_seed (liana6, swap_null_b, swap_b2, swap_t2) | 2 of 2 runs identical |
| AMENDMENT_4 sensitivity_permutations (liana6, swap_null_b, swap_b2, swap_t2) | 2 of 2 runs identical |
| AMENDMENT_4 sensitivity_permutations_calibration (liana6, swap_null_b, swap_b2, swap_t2) | 4 of 4 runs identical |

