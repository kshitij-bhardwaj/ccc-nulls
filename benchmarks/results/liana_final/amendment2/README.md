# AMENDMENT_2: robustness, held-out confirmation and reproducibility

Pre-registered in `benchmarks/liana_final/AMENDMENT_2.md` (committed before its runs). Run: `python -m benchmarks.liana_final.amendment2 --part all`. Measures as in AMENDMENT_1: false share in the top 1% (lower is better); power in the top 1% / AUROC planted vs one-sided / precision at k. Means over replicates [95% t-interval].

## A2.1 Held-out confirmation: Tabula Sapiens Lung

Calibration:

| | TS Lung, ligands shuffled (n=3) | TS Lung, receptors shuffled (n=3) |
|---|---|---|
| all rows | 0.966 | 0.939 |
| **liana6 (baseline)** | 0.833 [0.829, 0.837] | 0.880 [0.877, 0.883] |
| **swap: Null B** | 0.679 [0.668, 0.691] | 0.677 [0.666, 0.687] |
| **7th/8th: Null B** | 0.733 [0.725, 0.741] | 0.753 [0.741, 0.764] |
| swap: B2 | 0.738 [0.733, 0.744] | 0.758 [0.747, 0.769] |
| filter: detected ≥ 20% | 0.841 [0.838, 0.844] | 0.870 [0.867, 0.872] |
| filter: detected ≥ 30% | 0.855 [0.852, 0.859] | 0.865 [0.862, 0.869] |
| filter: both sides DE | 0.360 [0.352, 0.367] | 0.282 [0.281, 0.282] |
| liana6, mean aggregation | 0.820 [0.818, 0.821] | 0.859 [0.854, 0.864] |
| swap: Null B, mean aggregation | 0.690 [0.683, 0.698] | 0.665 [0.654, 0.676] |

Power (single pairs): power in the top 1% / AUROC / precision at k:

| | TS Lung, f = 1.5 (n=5) | TS Lung, f = 2.0 (n=5) |
|---|---|---|
| **liana6 (baseline)** | 0.99 / 0.945 / 0.14 | 1.00 / 0.945 / 0.11 |
| **swap: Null B** | 0.99 / 0.987 / 0.56 | 1.00 / 0.997 / 0.65 |
| **7th/8th: Null B** | 0.99 / 0.986 / 0.44 | 1.00 / 0.994 / 0.44 |
| swap: B2 | 0.92 / 0.951 / 0.39 | 0.99 / 0.972 / 0.47 |
| filter: detected ≥ 20% | 1.00 / 0.947 / 0.14 | 1.00 / 0.945 / 0.11 |
| filter: detected ≥ 30% | 0.92 / 0.928 / 0.19 | 1.00 / 0.956 / 0.11 |
| filter: both sides DE | 0.99 / 0.961 / 0.51 | 1.00 / 0.988 / 0.81 |
| liana6, mean aggregation | 0.99 / 0.955 / 0.32 | 1.00 / 0.974 / 0.37 |
| swap: Null B, mean aggregation | 0.99 / 0.989 / 0.58 | 1.00 / 0.998 / 0.69 |

Frozen rule (AMENDMENT_1, criteria 1 and 2):

| Variant | Calibration | Power | Verdict |
|---|---|---|---|
| **swap: Null B** | yes | yes | confirmed |
| **7th/8th: Null B** | yes | yes | confirmed |

## A2.2 Harder power (TS Mammary background)

Planted rows present: f = 1.1: complex 9.0, single 29.0, unequal 7.6; f = 1.25: complex 10.0, single 30.0, unequal 10.0; f = 1.5: complex 10.0, single 30.0, unequal 10.0.

**single** (30 pairs a replicate): power in the top 1% / AUROC planted vs one-sided / precision at k (all true rows):

| | TS Mammary, f = 1.1 (n=5) | TS Mammary, f = 1.25 (n=5) | TS Mammary, f = 1.5 (n=5) |
|---|---|---|---|
| **liana6 (baseline)** | 0.55 / 0.868 / 0.11 | 0.73 / 0.915 / 0.19 | 0.98 / 0.951 / 0.15 |
| **swap: Null B** | 0.65 / 0.914 / 0.21 | 0.99 / 0.985 / 0.54 | 1.00 / 0.996 / 0.72 |
| **7th/8th: Null B** | 0.59 / 0.898 / 0.14 | 0.89 / 0.960 / 0.44 | 0.99 / 0.992 / 0.52 |
| swap: B2 | 0.54 / 0.874 / 0.17 | 0.83 / 0.944 / 0.41 | 0.93 / 0.972 / 0.56 |
| filter: detected ≥ 20% | 0.28 / 0.666 / 0.08 | 0.77 / 0.906 / 0.19 | 0.98 / 0.952 / 0.15 |
| filter: detected ≥ 30% | 0.12 / 0.742 / 0.07 | 0.37 / 0.776 / 0.22 | 0.96 / 0.949 / 0.18 |
| filter: both sides DE | 0.57 / 0.874 / 0.16 | 0.86 / 0.950 / 0.49 | 0.99 / 0.987 / 0.81 |
| liana6, mean aggregation | 0.17 / 0.793 / 0.09 | 0.62 / 0.833 / 0.34 | 0.99 / 0.983 / 0.57 |
| swap: Null B, mean aggregation | 0.17 / 0.858 / 0.11 | 0.63 / 0.903 / 0.48 | 0.99 / 0.996 / 0.76 |

**complex** (10 pairs a replicate): power in the top 1% / AUROC planted vs one-sided / precision at k (all true rows):

| | TS Mammary, f = 1.1 (n=5) | TS Mammary, f = 1.25 (n=5) | TS Mammary, f = 1.5 (n=5) |
|---|---|---|---|
| **liana6 (baseline)** | 0.56 / 0.864 / 0.11 | 0.70 / 0.879 / 0.19 | 0.96 / 0.946 / 0.15 |
| **swap: Null B** | 0.62 / 0.912 / 0.21 | 0.98 / 0.974 / 0.54 | 0.98 / 0.996 / 0.72 |
| **7th/8th: Null B** | 0.60 / 0.898 / 0.14 | 0.88 / 0.937 / 0.44 | 1.00 / 0.992 / 0.52 |
| swap: B2 | 0.48 / 0.876 / 0.17 | 0.82 / 0.929 / 0.41 | 0.88 / 0.961 / 0.56 |
| filter: detected ≥ 20% | 0.10 / 0.598 / 0.08 | 0.68 / 0.849 / 0.19 | 0.96 / 0.947 / 0.15 |
| filter: detected ≥ 30% | 0.06 / 0.722 / 0.07 | 0.22 / 0.709 / 0.22 | 0.92 / 0.929 / 0.18 |
| filter: both sides DE | 0.58 / 0.874 / 0.16 | 0.86 / 0.930 / 0.49 | 0.98 / 0.985 / 0.81 |
| liana6, mean aggregation | 0.08 / 0.773 / 0.09 | 0.50 / 0.781 / 0.34 | 0.98 / 0.988 / 0.57 |
| swap: Null B, mean aggregation | 0.12 / 0.840 / 0.11 | 0.50 / 0.863 / 0.48 | 1.00 / 0.998 / 0.76 |

**unequal** (10 pairs a replicate): power in the top 1% / AUROC planted vs one-sided / precision at k (all true rows):

| | TS Mammary, f = 1.1 (n=5) | TS Mammary, f = 1.25 (n=5) | TS Mammary, f = 1.5 (n=5) |
|---|---|---|---|
| **liana6 (baseline)** | 0.30 / 0.722 / 0.11 | 0.38 / 0.790 / 0.19 | 0.62 / 0.838 / 0.15 |
| **swap: Null B** | 0.36 / 0.779 / 0.21 | 0.86 / 0.932 / 0.54 | 0.98 / 0.976 / 0.72 |
| **7th/8th: Null B** | 0.36 / 0.758 / 0.14 | 0.56 / 0.881 / 0.44 | 0.80 / 0.946 / 0.52 |
| swap: B2 | 0.22 / 0.774 / 0.17 | 0.62 / 0.868 / 0.41 | 0.82 / 0.922 / 0.56 |
| filter: detected ≥ 20% | 0.12 / 0.601 / 0.08 | 0.34 / 0.662 / 0.19 | 0.64 / 0.849 / 0.15 |
| filter: detected ≥ 30% | 0.06 / 0.636 / 0.07 | 0.18 / 0.659 / 0.22 | 0.48 / 0.740 / 0.18 |
| filter: both sides DE | 0.30 / 0.724 / 0.16 | 0.58 / 0.848 / 0.49 | 0.92 / 0.959 / 0.81 |
| liana6, mean aggregation | 0.10 / 0.672 / 0.09 | 0.26 / 0.677 / 0.34 | 0.62 / 0.809 / 0.57 |
| swap: Null B, mean aggregation | 0.10 / 0.757 / 0.11 | 0.34 / 0.817 / 0.48 | 0.72 / 0.904 / 0.76 |

Power gap to liana6 (single pairs, power in the top 1%):

| | f = 1.1 | f = 1.25 | f = 1.5 |
|---|---|---|---|
| **swap: Null B** | 0.10 [0.05, 0.15] | 0.26 [0.17, 0.35] | 0.02 [-0.02, 0.06] |
| **7th/8th: Null B** | 0.04 [-0.03, 0.11] | 0.16 [0.06, 0.26] | 0.01 [-0.02, 0.05] |
| filter: both sides DE | 0.01 [-0.02, 0.05] | 0.13 [0.07, 0.19] | 0.01 [-0.02, 0.05] |

## A2.3 Donor-aware permutations

(a) Labels shuffled within donor: share of rows (s ≠ t, score > 0) with p ≤ 0.05 (nominal 0.05), mean over replicates.

| Dataset (donors) | Tool | Null A, global | Null A, within donor | Null B, global | Null B, within donor |
|---|---|---|---|---|---|
| TS Lung (4) | CellPhoneDB | 0.236 [0.222, 0.250] | 0.064 [0.051, 0.077] | 0.060 [0.051, 0.070] | 0.004 [0.002, 0.005] |
| TS Lung (4) | CellChat | 0.297 [0.273, 0.321] | 0.073 [0.053, 0.093] | 0.075 [0.060, 0.091] | 0.004 [0.003, 0.005] |
| TNBC (10) | CellPhoneDB | 0.369 [0.364, 0.374] | 0.062 [0.047, 0.076] | 0.128 [0.120, 0.136] | 0.003 [0.002, 0.005] |
| TNBC (10) | CellChat | 0.506 [0.487, 0.524] | 0.060 [0.041, 0.078] | 0.184 [0.173, 0.196] | 0.003 [0.002, 0.004] |
| TS Mammary (2) | CellPhoneDB | 0.178 [0.141, 0.216] | 0.075 [0.056, 0.094] | 0.037 [0.028, 0.046] | 0.004 [0.003, 0.005] |
| TS Mammary (2) | CellChat | 0.193 [0.154, 0.232] | 0.102 [0.073, 0.130] | 0.040 [0.034, 0.046] | 0.006 [0.004, 0.007] |
| TS Liver (4) | CellPhoneDB | 0.304 [0.279, 0.330] | 0.060 [0.047, 0.073] | 0.107 [0.097, 0.117] | 0.003 [0.002, 0.003] |
| TS Liver (4) | CellChat | 0.377 [0.332, 0.422] | 0.075 [0.051, 0.100] | 0.146 [0.124, 0.169] | 0.004 [0.002, 0.005] |

(b) Real labels: swap_null_b with within-donor against global permutations.

| Dataset | Top-1% Jaccard | Calls, global | Calls, within donor | Call Jaccard |
|---|---|---|---|---|
| TNBC | 0.746 | 24,005 | 23,393 | 0.913 |
| TS Mammary | 0.799 | 13,208 | 12,967 | 0.913 |
| TS Liver | 0.597 | 33,740 | 33,153 | 0.868 |

## A2.4 Null B against simple filters

Settings (calibration scenarios and power settings of this amendment) in which the Null B variant beats the filter: lower false share in the top 1% / higher AUROC (single pairs).

| | filter: detected ≥ 20% | filter: detected ≥ 30% | filter: both sides DE |
|---|---|---|---|
| **swap: Null B** | calibration 18/18; AUROC 9/9 | calibration 18/18; AUROC 9/9 | calibration 0/18; AUROC 9/9 |
| **7th/8th: Null B** | calibration 18/18; AUROC 9/9 | calibration 18/18; AUROC 9/9 | calibration 0/18; AUROC 7/9 |

## A2.5 Sensitivity (TS Mammary)

Seed (real labels, against seed 1337): Spearman of specificity_rank / Jaccard of the top 1%.

| Seed | liana6 | swap: Null B | 7th/8th: Null B |
|---|---|---|---|
| 1 | 0.9979 / 0.964 | 0.9965 / 0.935 | 0.9947 / 0.961 |
| 2 | 0.9977 / 0.979 | 0.9960 / 0.973 | 0.9943 / 0.964 |

Permutations, 1,000 vs 5,000 (the hypothalamus; see Deviations):

| Real labels | liana6 | swap: Null B | 7th/8th: Null B |
|---|---|---|---|
| Spearman / top-1% Jaccard | 0.9991 / 1.000 | 0.9969 / 1.000 | 0.9975 / 0.889 |

| | Hypothalamus @n_perms=1000, receptors shuffled (n=1) | Hypothalamus @n_perms=5000, receptors shuffled (n=1) | TS Mammary @n_perms=1000, receptors shuffled (n=1) | TS Mammary @n_perms=5000, receptors shuffled (n=1) |
|---|---|---|---|---|
| all rows | 0.882 | 0.882 | 0.919 | 0.919 |
| **liana6 (baseline)** | 0.758 | 0.758 | 0.804 | 0.806 |
| **swap: Null B** | 0.192 | 0.183 | 0.694 | 0.695 |
| **7th/8th: Null B** | 0.208 | 0.225 | 0.673 | 0.653 |

Calibration under each setting (rep 0):

| | TS Mammary @cells=0.25, receptors shuffled (n=1) | TS Mammary @cells=0.5, receptors shuffled (n=1) | TS Mammary @expr_prop=0.05, ligands shuffled (n=1) | TS Mammary @expr_prop=0.05, receptors shuffled (n=1) | TS Mammary @expr_prop=0.2, ligands shuffled (n=1) | TS Mammary @expr_prop=0.2, receptors shuffled (n=1) | TS Mammary @labels=broad_cell_class, ligands shuffled (n=1) | TS Mammary @labels=broad_cell_class, receptors shuffled (n=1) | TS Mammary @labels=free_annotation, ligands shuffled (n=1) | TS Mammary @labels=free_annotation, receptors shuffled (n=1) |
|---|---|---|---|---|---|---|---|---|---|---|
| all rows | 0.921 | 0.919 | 0.974 | 0.923 | 0.975 | 0.919 | 0.974 | 0.917 | 0.974 | 0.924 |
| **liana6 (baseline)** | 0.771 | 0.815 | 0.874 | 0.849 | 0.752 | 0.734 | 0.853 | 0.749 | 0.839 | 0.802 |
| **swap: Null B** | 0.655 | 0.688 | 0.768 | 0.744 | 0.509 | 0.566 | 0.719 | 0.632 | 0.684 | 0.694 |
| **7th/8th: Null B** | 0.713 | 0.700 | 0.812 | 0.726 | 0.555 | 0.535 | 0.748 | 0.642 | 0.722 | 0.736 |
| swap: B2 | 0.695 | 0.730 | 0.844 | 0.789 | 0.649 | 0.549 | 0.809 | 0.656 | 0.763 | 0.731 |
| filter: detected ≥ 20% | 0.783 | 0.826 | 0.879 | 0.843 | 0.752 | 0.732 | 0.881 | 0.768 | 0.847 | 0.803 |
| filter: detected ≥ 30% | 0.801 | 0.857 | 0.889 | 0.849 | 0.799 | 0.746 | 0.897 | 0.763 | 0.866 | 0.812 |
| filter: both sides DE | 0.387 | 0.451 | 0.682 | 0.530 | 0.356 | 0.064 | 0.506 | 0.201 | 0.554 | 0.424 |
| liana6, mean aggregation | 0.778 | 0.773 | 0.851 | 0.798 | 0.728 | 0.726 | 0.849 | 0.740 | 0.809 | 0.786 |
| swap: Null B, mean aggregation | 0.659 | 0.689 | 0.819 | 0.749 | 0.519 | 0.542 | 0.778 | 0.635 | 0.717 | 0.703 |

Power under each setting (f = 1.5, single pairs):

| | TS Mammary @cells=0.25, f = 1.5 (n=1) | TS Mammary @cells=0.5, f = 1.5 (n=1) | TS Mammary @expr_prop=0.05, f = 1.5 (n=2) | TS Mammary @expr_prop=0.2, f = 1.5 (n=2) |
|---|---|---|---|---|
| **liana6 (baseline)** | 0.90 / 0.920 / 0.12 | 0.97 / 0.949 / 0.12 | 1.00 / 0.958 / 0.17 | 0.95 / 0.924 / 0.17 |
| **swap: Null B** | 0.93 / 0.977 / 0.64 | 1.00 / 0.993 / 0.60 | 1.00 / 0.997 / 0.70 | 1.00 / 0.994 / 0.70 |
| **7th/8th: Null B** | 0.93 / 0.974 / 0.50 | 1.00 / 0.988 / 0.45 | 1.00 / 0.994 / 0.41 | 1.00 / 0.990 / 0.48 |
| swap: B2 | 0.93 / 0.946 / 0.48 | 1.00 / 0.970 / 0.50 | 0.98 / 0.969 / 0.46 | 0.88 / 0.955 / 0.48 |
| filter: detected ≥ 20% | 0.90 / 0.925 / 0.12 | 0.97 / 0.949 / 0.12 | 1.00 / 0.958 / 0.17 | 0.95 / 0.924 / 0.17 |
| filter: detected ≥ 30% | 0.90 / 0.913 / 0.14 | 0.93 / 0.928 / 0.12 | 1.00 / 0.965 / 0.19 | 0.95 / 0.935 / 0.18 |
| filter: both sides DE | 0.90 / 0.931 / 0.29 | 0.97 / 0.972 / 0.50 | 1.00 / 0.995 / 0.89 | 1.00 / 0.990 / 0.88 |
| liana6, mean aggregation | 0.90 / 0.927 / 0.55 | 1.00 / 0.982 / 0.62 | 1.00 / 0.984 / 0.59 | 0.97 / 0.974 / 0.59 |
| swap: Null B, mean aggregation | 0.93 / 0.970 / 0.67 | 1.00 / 0.996 / 0.67 | 1.00 / 0.998 / 0.77 | 1.00 / 0.998 / 0.77 |

## A2.6 The p-values themselves

Calibration rerun with the same seeds (Mammary, Liver):

| | TS Liver, ligands shuffled (n=3) | TS Liver, receptors shuffled (n=3) | TS Mammary, ligands shuffled (n=3) | TS Mammary, receptors shuffled (n=3) |
|---|---|---|---|---|
| all rows | 0.963 | 0.929 | 0.974 | 0.920 |
| **liana6 (baseline)** | 0.850 [0.841, 0.858] | 0.827 [0.817, 0.836] | 0.857 [0.853, 0.862] | 0.799 [0.782, 0.815] |
| **swap: Null B** | 0.628 [0.590, 0.666] | 0.602 [0.557, 0.646] | 0.677 [0.668, 0.686] | 0.694 [0.669, 0.718] |
| **7th/8th: Null B** | 0.693 [0.685, 0.701] | 0.690 [0.673, 0.708] | 0.726 [0.710, 0.742] | 0.696 [0.624, 0.769] |
| swap: B2 | 0.712 [0.663, 0.761] | 0.692 [0.657, 0.727] | 0.779 [0.765, 0.792] | 0.733 [0.721, 0.746] |
| filter: detected ≥ 20% | 0.858 [0.849, 0.866] | 0.823 [0.815, 0.832] | 0.855 [0.850, 0.860] | 0.812 [0.801, 0.824] |
| filter: detected ≥ 30% | 0.903 [0.895, 0.911] | 0.827 [0.814, 0.839] | 0.871 [0.869, 0.873] | 0.823 [0.809, 0.836] |
| filter: both sides DE | 0.376 [0.357, 0.395] | 0.195 [0.186, 0.205] | 0.592 [0.586, 0.599] | 0.406 [0.389, 0.423] |
| liana6, mean aggregation | 0.806 [0.798, 0.814] | 0.822 [0.795, 0.848] | 0.802 [0.796, 0.807] | 0.771 [0.745, 0.797] |
| swap: Null B, mean aggregation | 0.659 [0.636, 0.682] | 0.577 [0.529, 0.626] | 0.719 [0.715, 0.723] | 0.693 [0.647, 0.738] |

Reproducibility: 12 of 12 reruns give the false shares of `calibration.json` exactly.

Share of false rows (shuffled side) with p ≤ α, pooled over runs; BH at 5% within each tool: share of discoveries that are false rows (beside the share of false rows overall).

| Tool | Null | p ≤ 0.01 | p ≤ 0.05 | p ≤ 0.1 | BH discoveries | false share of discoveries | false share of all rows |
|---|---|---|---|---|---|---|---|
| CellPhoneDB | Null A | 0.3369 | 0.3833 | 0.4135 | 70,597 | 0.939 | 0.946 |
| CellPhoneDB | Null B | 0.0055 | 0.0262 | 0.0524 | 0 | nan | 0.946 |
| CellPhoneDB | B2 | 0.0000 | 0.0009 | 0.0029 | 2 | 0.000 | 0.946 |
| CellChat | Null A | 0.1689 | 0.1876 | 0.1992 | 34,106 | 0.927 | 0.946 |
| CellChat | Null B | 0.0026 | 0.0130 | 0.0264 | 0 | nan | 0.946 |
| CellChat | B2 | 0.0000 | 0.0009 | 0.0029 | 2 | 0.000 | 0.946 |

## A2.7 A second species: Allen brain (mouse)

| | Allen brain, ligands shuffled (n=3) | Allen brain, receptors shuffled (n=3) |
|---|---|---|
| all rows | 0.930 | 0.913 |
| **liana6 (baseline)** | 0.925 [0.925, 0.926] | 0.836 [0.831, 0.840] |
| **swap: Null B** | 0.825 [0.821, 0.829] | 0.726 [0.708, 0.744] |
| **7th/8th: Null B** | 0.808 [0.776, 0.839] | 0.690 [0.656, 0.723] |
| swap: B2 | 0.853 [0.849, 0.857] | 0.735 [0.725, 0.745] |
| filter: detected ≥ 20% | 0.915 [0.913, 0.917] | 0.827 [0.823, 0.831] |
| filter: detected ≥ 30% | 0.890 [0.889, 0.891] | 0.837 [0.834, 0.841] |
| filter: both sides DE | 0.460 [0.457, 0.463] | 0.391 [0.388, 0.394] |
| liana6, mean aggregation | 0.924 [0.923, 0.926] | 0.802 [0.795, 0.809] |
| swap: Null B, mean aggregation | 0.806 [0.803, 0.809] | 0.710 [0.674, 0.746] |

AMENDMENT_1's criterion 1 here: **swap: Null B** yes, **7th/8th: Null B** yes

## A2.8 Biology: every test, multiplicity and detectable effects

Odds ratio among the baseline's calls, kept vs dropped [95% cluster-bootstrap interval]; two-sided bootstrap p; BH q within the variant and across all; detectable: the odds ratios each test could detect with 80% power (α = 0.05).

****swap: Null B****

| Test | OR [95% CI] | p | q (variant) | q (all) | detectable |
|---|---|---|---|---|---|
| tnbc_wu2021 spatial | 0.50 [0.29, 1.44] | 0.1540 | 0.308 | 0.366 | ≤ 0.33 or ≥ 2.99 |
| allen_brain_atlas spatial | 1.35 [0.96, 1.86] | 0.0920 | 0.276 | 0.291 | ≤ 0.62 or ≥ 1.61 |
| moffitt2018_hypothalamus spatial | 1.15 [0.72, 2.13] | 0.5950 | 0.595 | 0.797 | ≤ 0.46 or ≥ 2.18 |
| tnbc_wu2021 cytosig | 0.51 [0.21, 1.40] | 0.1530 | 0.308 | 0.366 | ≤ 0.27 or ≥ 3.75 |
| tnbc_wu2021 progeny | 6.28 [2.24, 29.91] | 0.0020 | 0.012 | 0.012 | ≤ 0.18 or ≥ 5.62 |
| allen_brain_atlas cytosig | 5.53 [2.80, 15.93] | 0.0005 | 0.006 | 0.004 | ≤ 0.30 or ≥ 3.36 |
| allen_brain_atlas progeny | 0.63 [0.19, 5.41] | 0.5260 | 0.595 | 0.758 | ≤ 0.10 or ≥ 9.65 |
| tabula_sapiens_mammary cytosig | 1.51 [0.63, 16.75] | 0.3430 | 0.588 | 0.602 | ≤ 0.10 or ≥ 9.60 |
| tabula_sapiens_mammary progeny | 0.77 [0.47, 1.52] | 0.4530 | 0.595 | 0.742 | ≤ 0.42 or ≥ 2.40 |
| tabula_sapiens_liver cytosig | 0.79 [0.32, 2.64] | 0.5840 | 0.595 | 0.797 | ≤ 0.25 or ≥ 4.03 |
| tabula_sapiens_liver progeny | 0.47 [0.21, 0.95] | 0.0340 | 0.136 | 0.128 | ≤ 0.35 or ≥ 2.86 |
| moffitt2018_hypothalamus cytosig | – (keeps or drops none) | | | | |
| moffitt2018_hypothalamus progeny | 0.66 [0.02, 2.90] | 0.4790 | 0.595 | 0.742 | ≤ 0.02 or ≥ 57.42 |

****7th/8th: Null B****

| Test | OR [95% CI] | p | q (variant) | q (all) | detectable |
|---|---|---|---|---|---|
| tnbc_wu2021 spatial | 1.01 [0.48, 1.97] | 0.9960 | 0.996 | 0.996 | ≤ 0.36 or ≥ 2.77 |
| allen_brain_atlas spatial | 1.07 [0.81, 1.35] | 0.6620 | 0.717 | 0.813 | ≤ 0.69 or ≥ 1.46 |
| moffitt2018_hypothalamus spatial | 1.10 [0.72, 1.88] | 0.6450 | 0.717 | 0.813 | ≤ 0.51 or ≥ 1.95 |
| tnbc_wu2021 cytosig | 0.59 [0.23, 1.86] | 0.3060 | 0.663 | 0.549 | ≤ 0.23 or ≥ 4.31 |
| tnbc_wu2021 progeny | 3.77 [1.15, 15.14] | 0.0330 | 0.107 | 0.128 | ≤ 0.16 or ≥ 6.40 |
| allen_brain_atlas cytosig | 1.83 [1.15, 3.34] | 0.0160 | 0.082 | 0.079 | ≤ 0.45 or ≥ 2.21 |
| allen_brain_atlas progeny | 0.64 [0.19, 5.63] | 0.5140 | 0.717 | 0.758 | ≤ 0.11 or ≥ 9.50 |
| tabula_sapiens_mammary cytosig | 0.76 [0.35, 4.16] | 0.6580 | 0.717 | 0.813 | ≤ 0.18 or ≥ 5.50 |
| tabula_sapiens_mammary progeny | 0.60 [0.35, 1.41] | 0.1670 | 0.434 | 0.366 | ≤ 0.37 or ≥ 2.67 |
| tabula_sapiens_liver cytosig | 0.78 [0.27, 3.01] | 0.6260 | 0.717 | 0.811 | ≤ 0.19 or ≥ 5.16 |
| tabula_sapiens_liver progeny | 0.24 [0.09, 0.70] | 0.0190 | 0.082 | 0.083 | ≤ 0.23 or ≥ 4.39 |
| moffitt2018_hypothalamus cytosig | 0.02 [0.01, 0.23] | 0.0005 | 0.007 | 0.004 | ≤ 0.08 or ≥ 12.71 |
| moffitt2018_hypothalamus progeny | 0.66 [0.02, 2.90] | 0.4790 | 0.717 | 0.742 | ≤ 0.02 or ≥ 57.42 |

**swap: Null B, donor permutations**

| Test | OR [95% CI] | p | q (variant) | q (all) | detectable |
|---|---|---|---|---|---|
| tnbc_wu2021 spatial | 0.60 [0.39, 1.81] | 0.2560 | 0.256 | 0.493 | ≤ 0.34 or ≥ 2.97 |
| tnbc_wu2021 cytosig | 0.50 [0.19, 1.49] | 0.1780 | 0.256 | 0.380 | ≤ 0.24 or ≥ 4.18 |
| tnbc_wu2021 progeny | 6.17 [2.38, 18.36] | 0.0010 | 0.003 | 0.007 | ≤ 0.23 or ≥ 4.35 |

**swap: B2**

| Test | OR [95% CI] | p | q (variant) | q (all) | detectable |
|---|---|---|---|---|---|
| tnbc_wu2021 spatial | 1.33 [0.80, 2.22] | 0.2200 | 0.318 | 0.434 | ≤ 0.49 or ≥ 2.04 |
| allen_brain_atlas spatial | 1.48 [0.96, 2.39] | 0.0710 | 0.208 | 0.255 | ≤ 0.52 or ≥ 1.92 |
| moffitt2018_hypothalamus spatial | 1.82 [0.93, 4.33] | 0.0800 | 0.208 | 0.263 | ≤ 0.33 or ≥ 3.04 |
| tnbc_wu2021 cytosig | 0.44 [0.18, 2.44] | 0.2020 | 0.318 | 0.411 | ≤ 0.17 or ≥ 5.95 |
| tnbc_wu2021 progeny | 4.43 [1.21, 29.01] | 0.0190 | 0.082 | 0.083 | ≤ 0.09 or ≥ 10.81 |
| allen_brain_atlas cytosig | 2.94 [1.50, 7.04] | 0.0020 | 0.013 | 0.012 | ≤ 0.34 or ≥ 2.92 |
| allen_brain_atlas progeny | 0.47 [0.12, 3.44] | 0.2720 | 0.354 | 0.500 | ≤ 0.11 or ≥ 8.82 |
| tabula_sapiens_mammary cytosig | 1.16 [0.48, 10.54] | 0.7000 | 0.700 | 0.830 | ≤ 0.11 or ≥ 9.48 |
| tabula_sapiens_mammary progeny | 0.71 [0.35, 2.91] | 0.6690 | 0.700 | 0.813 | ≤ 0.21 or ≥ 4.66 |
| tabula_sapiens_liver cytosig | 0.61 [0.28, 7.92] | 0.3970 | 0.469 | 0.667 | ≤ 0.13 or ≥ 7.75 |
| tabula_sapiens_liver progeny | 0.42 [0.15, 2.66] | 0.2030 | 0.318 | 0.411 | ≤ 0.15 or ≥ 6.72 |
| moffitt2018_hypothalamus cytosig | 0.16 [0.01, 0.51] | 0.0010 | 0.013 | 0.007 | ≤ 0.05 or ≥ 20.48 |
| moffitt2018_hypothalamus progeny | 0.38 [0.02, 1.51] | 0.1600 | 0.318 | 0.366 | ≤ 0.03 or ≥ 37.84 |

**7th/8th: B2**

| Test | OR [95% CI] | p | q (variant) | q (all) | detectable |
|---|---|---|---|---|---|
| tnbc_wu2021 spatial | 1.15 [0.62, 1.95] | 0.6250 | 0.981 | 0.811 | ≤ 0.44 or ≥ 2.25 |
| allen_brain_atlas spatial | 1.24 [0.93, 1.64] | 0.1480 | 0.347 | 0.366 | ≤ 0.66 or ≥ 1.51 |
| moffitt2018_hypothalamus spatial | 1.36 [0.78, 2.81] | 0.2660 | 0.494 | 0.500 | ≤ 0.40 or ≥ 2.50 |
| tnbc_wu2021 cytosig | 0.95 [0.33, 5.87] | 0.9190 | 0.981 | 0.968 | ≤ 0.13 or ≥ 7.43 |
| tnbc_wu2021 progeny | 3.96 [1.21, 25.44] | 0.0260 | 0.113 | 0.108 | ≤ 0.12 or ≥ 8.40 |
| allen_brain_atlas cytosig | 1.78 [1.13, 3.50] | 0.0090 | 0.058 | 0.051 | ≤ 0.45 or ≥ 2.24 |
| allen_brain_atlas progeny | 0.98 [0.42, 5.82] | 0.9810 | 0.981 | 0.994 | ≤ 0.15 or ≥ 6.81 |
| tabula_sapiens_mammary cytosig | 0.71 [0.29, 6.32] | 0.7040 | 0.981 | 0.830 | ≤ 0.12 or ≥ 8.58 |
| tabula_sapiens_mammary progeny | 0.39 [0.19, 1.92] | 0.1520 | 0.347 | 0.366 | ≤ 0.20 or ≥ 4.94 |
| tabula_sapiens_liver cytosig | 1.14 [0.30, 11.30] | 0.8830 | 0.981 | 0.960 | ≤ 0.09 or ≥ 11.27 |
| tabula_sapiens_liver progeny | 0.74 [0.25, 4.39] | 0.7590 | 0.981 | 0.859 | ≤ 0.12 or ≥ 8.10 |
| moffitt2018_hypothalamus cytosig | 0.02 [0.00, 0.19] | 0.0005 | 0.007 | 0.004 | ≤ 0.08 or ≥ 13.07 |
| moffitt2018_hypothalamus progeny | 0.38 [0.02, 1.51] | 0.1600 | 0.347 | 0.366 | ≤ 0.03 or ≥ 37.84 |

**7th/8th: Null A**

| Test | OR [95% CI] | p | q (variant) | q (all) | detectable |
|---|---|---|---|---|---|
| tnbc_wu2021 spatial | 0.95 [0.34, 6.02] | 0.9630 | 0.963 | 0.988 | ≤ 0.12 or ≥ 8.43 |
| allen_brain_atlas spatial | 0.51 [0.36, 0.85] | 0.0140 | 0.056 | 0.074 | ≤ 0.53 or ≥ 1.88 |
| moffitt2018_hypothalamus spatial | 0.88 [0.44, 3.03] | 0.7610 | 0.963 | 0.859 | ≤ 0.26 or ≥ 3.91 |
| tnbc_wu2021 cytosig | 0.90 [0.27, 3.73] | 0.9100 | 0.963 | 0.968 | ≤ 0.16 or ≥ 6.33 |
| tnbc_wu2021 progeny | 1.31 [0.26, 12.15] | 0.7380 | 0.963 | 0.857 | ≤ 0.06 or ≥ 16.33 |
| allen_brain_atlas cytosig | 0.26 [0.15, 0.50] | 0.0005 | 0.003 | 0.004 | ≤ 0.44 or ≥ 2.30 |
| allen_brain_atlas progeny | 0.97 [0.20, 3.69] | 0.9470 | 0.963 | 0.984 | ≤ 0.14 or ≥ 7.19 |
| tabula_sapiens_mammary cytosig | 0.27 [0.04, 1.13] | 0.0780 | 0.234 | 0.263 | ≤ 0.08 or ≥ 11.88 |
| tabula_sapiens_mammary progeny | 1.88 [0.25, 5.46] | 0.5450 | 0.963 | 0.769 | ≤ 0.11 or ≥ 9.20 |
| tabula_sapiens_liver cytosig | 1.13 [0.30, 7.65] | 0.8870 | 0.963 | 0.960 | ≤ 0.10 or ≥ 9.54 |
| tabula_sapiens_liver progeny | 0.39 [0.11, 3.47] | 0.3540 | 0.850 | 0.608 | ≤ 0.09 or ≥ 11.18 |
| moffitt2018_hypothalamus cytosig | 0.02 [0.01, 0.23] | 0.0005 | 0.003 | 0.004 | ≤ 0.08 or ≥ 12.71 |
| moffitt2018_hypothalamus progeny | – (keeps or drops none) | | | | |

**LIANA+ as shipped**

| Test | OR [95% CI] | p | q (variant) | q (all) | detectable |
|---|---|---|---|---|---|
| tnbc_wu2021 spatial | 0.88 [0.57, 1.53] | 0.5280 | 0.624 | 0.758 | ≤ 0.50 or ≥ 1.99 |
| allen_brain_atlas spatial | 1.58 [0.86, 3.47] | 0.1440 | 0.271 | 0.366 | ≤ 0.37 or ≥ 2.74 |
| moffitt2018_hypothalamus spatial | 1.31 [0.51, 2.66] | 0.5030 | 0.624 | 0.758 | ≤ 0.31 or ≥ 3.19 |
| tnbc_wu2021 cytosig | 1.11 [0.72, 2.28] | 0.5860 | 0.635 | 0.797 | ≤ 0.43 or ≥ 2.33 |
| tnbc_wu2021 progeny | 2.64 [1.64, 7.44] | 0.0005 | 0.002 | 0.004 | ≤ 0.34 or ≥ 2.92 |
| allen_brain_atlas cytosig | 5.52 [2.73, 33.93] | 0.0005 | 0.002 | 0.004 | ≤ 0.18 or ≥ 5.46 |
| allen_brain_atlas progeny | 1.23 [0.24, 6.00] | 0.8440 | 0.844 | 0.939 | ≤ 0.11 or ≥ 9.40 |
| tabula_sapiens_mammary cytosig | 8.34 [2.94, 38.16] | 0.0005 | 0.002 | 0.004 | ≤ 0.15 or ≥ 6.87 |
| tabula_sapiens_mammary progeny | 2.63 [0.68, 22.09] | 0.1560 | 0.271 | 0.366 | ≤ 0.10 or ≥ 10.23 |
| tabula_sapiens_liver cytosig | 1.29 [0.90, 3.57] | 0.1670 | 0.271 | 0.366 | ≤ 0.37 or ≥ 2.69 |
| tabula_sapiens_liver progeny | 10.82 [3.93, 27.06] | 0.0005 | 0.002 | 0.004 | ≤ 0.24 or ≥ 4.14 |
| moffitt2018_hypothalamus cytosig | 2.17 [0.18, 6.18] | 0.4730 | 0.624 | 0.742 | ≤ 0.08 or ≥ 13.22 |
| moffitt2018_hypothalamus progeny | 0.11 [0.03, 1.21] | 0.1580 | 0.271 | 0.366 | ≤ 0.04 or ≥ 27.58 |

## A2.9 What Null B drops, and why

Among liana6's calls (s ≠ t): rows swap_null_b drops vs keeps. Limiting side: the side with the larger Null B p. Breadth: share of cell types in which the limiting side's gene is detected in ≥ 10% of cells. Mann–Whitney: dropped broader than kept.

| Dataset | Tool | Dropped | Kept | Limiting side of dropped (ligand / receptor) | Breadth, dropped vs kept (median) | p |
|---|---|---|---|---|---|---|
| TNBC | CellPhoneDB | 3,926 | 18,356 | 1,814 / 2,101 | 0.90 vs 0.59 | 0 |
| TNBC | CellChat | 3,926 | 18,356 | 1,585 / 1,714 | 0.90 vs 0.55 | 0 |
| Allen brain | CellPhoneDB | 9,443 | 31,456 | 4,058 / 5,317 | 0.91 vs 0.83 | 8.6e-167 |
| Allen brain | CellChat | 9,443 | 31,456 | 3,962 / 5,278 | 0.91 vs 0.87 | 3e-132 |
| TS Mammary | CellPhoneDB | 3,468 | 11,178 | 1,709 / 1,743 | 0.88 vs 0.75 | 0.0031 |
| TS Mammary | CellChat | 3,468 | 11,178 | 1,197 / 1,191 | 0.88 vs 0.75 | 1.9e-15 |
| TS Liver | CellPhoneDB | 4,394 | 26,979 | 2,619 / 1,766 | 0.96 vs 0.92 | 2.9e-92 |
| TS Liver | CellChat | 4,394 | 26,979 | 2,586 / 1,761 | 0.96 vs 0.92 | 2.6e-85 |
| Hypothalamus | CellPhoneDB | 146 | 1,083 | 110 / 36 | 1.00 vs 0.75 | 1e-21 |
| Hypothalamus | CellChat | 146 | 1,083 | 96 / 30 | 1.00 vs 0.75 | 2e-23 |

Top dropped and promoted rows: `cases_<dataset>_<tool>_dropped.csv`, `cases_<dataset>_promoted.csv`.

## A2.10 Run time

| Dataset | Cells | li.mt.rank_aggregate (s) | This pipeline (s) | Ratio |
|---|---|---|---|---|
| TS Mammary | 30,934 | 399 | 2,208 | 5.54 |
| Hypothalamus | 29,917 | 201 | 1,562 | 7.78 |

This pipeline: six methods (LIANA+'s five and CellChat), Null B, and B2 with its decoy genes.

## Deviations from AMENDMENT_2

- **Permutations (A2.5).** 5,000 permutations on Tabula Sapiens Mammary need a 2 × 5,000 × 86,000 permutation tensor (about 7 GB), beyond the 8 GB machine. The comparison was run on the hypothalamus (10,197 rows), real labels and the receptor-shuffled calibration (rep 0).
- **Lung (A2.1, A2.3): another machine, not a subsample.** A full-Lung run was killed by the operating system for memory, alone on the 8 GB machine, and so was the pre-registered fallback of 30,000 cells: memory is set by the permutation tensor, 2 × 1,000 × rows, and rows grow with the square of the number of cell types (34 in Lung, about 370,000 rows), not with the number of cells. Instead of a subsample, every Lung run used all 65,847 cells on https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-donor-liver, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-donor-mammary, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-donor-tnbc, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-donor-0, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-donor-1, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-donor-2, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-global-0, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-global-1, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-global-2, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-ligand-0, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-ligand-1, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-ligand-2, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-power, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-receptor-0, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-receptor-1, https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-receptor-2, in the same locked environment (`requirements/lock-liana-final.txt`), with the same code, seeds and settings; its result files were merged with `benchmarks/liana_final/merge_results.py`, and each Lung run records the machine.

