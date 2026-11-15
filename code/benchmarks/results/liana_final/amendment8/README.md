# AMENDMENT_8 results: NATMI's specificity replaced too, through every test

Plan: [`AMENDMENT_8.md`](../../../liana_final/AMENDMENT_8.md). `swap_<n>`: CellPhoneDB's and CellChat's p replaced by null n's; `swap3_<n>`: Connectome's too; `swap4_<n>`: NATMI's too. Names as in the paper: JLP = Null A, BPT = Null B, MDT = B2, LOO-t = T2. Intervals: 95% t-intervals over replicates.

Where it ran: 31 of the 42 sessions in `kaggle/kernels.py` on Kaggle; the other 11 (`a8-sensitivity-*`, `a8-donor-mammary`, `a8-donor-tnbc`, `a8-tool-real-*`) on a lab server, from the same commit (be09388) and lock file, one CPU each; every run records its machine. Biology and case studies ran locally from the real-data tables. The tool-level spatial q-values are recomputed here over every data set (each session wrote q over its own).

## Calibration: false share of the top 1% (lower is better)

swap3: CellPhoneDB's, CellChat's and Connectome's p replaced; swap4: NATMI's too. Same runs, so paired.

| Setting | reps | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|---|
| Mammary, receptors | 3 | 0.799 [0.782, 0.815] | 0.499 [0.459, 0.539] | 0.393 [0.352, 0.434] | 0.694 [0.678, 0.710] | 0.677 [0.628, 0.727] | 0.832 [0.818, 0.847] | 0.873 [0.869, 0.878] | 0.845 [0.832, 0.858] |
| Mammary, ligands | 3 | 0.857 [0.853, 0.862] | 0.600 [0.579, 0.620] | 0.553 [0.546, 0.560] | 0.784 [0.774, 0.794] | 0.797 [0.772, 0.822] | 0.929 [0.924, 0.934] | 0.954 [0.944, 0.964] | 0.931 [0.930, 0.931] |
| Liver, receptors | 3 | 0.827 [0.817, 0.836] | 0.322 [0.300, 0.344] | 0.395 [0.338, 0.451] | 0.626 [0.568, 0.684] | 0.580 [0.544, 0.617] | 0.820 [0.788, 0.853] | 0.868 [0.848, 0.889] | 0.873 [0.868, 0.879] |
| Liver, ligands | 3 | 0.850 [0.841, 0.858] | 0.467 [0.439, 0.495] | 0.453 [0.378, 0.528] | 0.764 [0.730, 0.798] | 0.791 [0.770, 0.812] | 0.903 [0.897, 0.910] | 0.922 [0.915, 0.928] | 0.955 [0.953, 0.956] |
| Lung (held out), receptors | 3 | 0.880 [0.877, 0.883] | 0.434 [0.431, 0.438] | 0.448 [0.414, 0.482] | 0.701 [0.677, 0.725] | 0.709 [0.691, 0.726] | 0.875 [0.872, 0.878] | 0.925 [0.924, 0.926] | 0.928 [0.923, 0.933] |
| Lung (held out), ligands | 3 | 0.833 [0.829, 0.837] | 0.494 [0.485, 0.503] | 0.452 [0.305, 0.599] | 0.739 [0.730, 0.748] | 0.750 [0.740, 0.759] | 0.911 [0.910, 0.912] | 0.935 [0.933, 0.937] | 0.963 [0.959, 0.967] |
| Allen brain, receptors | 3 | 0.836 [0.831, 0.840] | 0.467 [0.459, 0.475] | 0.429 [0.422, 0.436] | 0.562 [0.547, 0.577] | 0.509 [0.459, 0.560] | 0.897 [0.893, 0.901] | 0.928 [0.916, 0.940] | 0.866 [0.863, 0.869] |
| Allen brain, ligands | 3 | 0.925 [0.925, 0.926] | 0.532 [0.523, 0.540] | 0.450 [0.447, 0.454] | 0.744 [0.716, 0.771] | 0.713 [0.675, 0.752] | 0.952 [0.949, 0.956] | 0.945 [0.943, 0.947] | 0.877 [0.872, 0.881] |

## AMENDMENT_1's rule for swap4 (and swap3, from the same runs)

Criterion 1: false share of the top 1% below liana6's on both sides; criterion 2: power in the top 1% ≥ 0.8 and AUROC ≥ liana6's at f = 2.0.

| Variant | 1, Mammary / Liver | 2, Mammary | 1, Lung | 2, Lung | Held-out Lung |
|---|---|---|---|---|---|
| swap4_null_a | fails / fails | fails | fails | fails | not confirmed (calibration, power) |
| swap4_null_b | holds / holds | holds | holds | holds | confirmed |
| swap4_b2 | holds / holds | fails | holds | fails | not confirmed (power) |
| swap4_t2 | fails / fails | holds | fails | holds | not confirmed (calibration) |
| swap3_null_a | fails / fails | fails | fails | fails | not confirmed (calibration, power) |
| swap3_null_b | holds / holds | holds | holds | holds | confirmed |
| swap3_b2 | holds / holds | holds | holds | fails | not confirmed (power) |
| swap3_t2 | fails / fails | holds | fails | holds | not confirmed (calibration) |

## Power, F2 design: power in the top 1% / AUROC, single pairs

| Effect | reps | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|---|
| f = 1.5 | 5 | 1.00 / 0.947 | 1.00 / 0.999 | 1.00 / 1.000 | 0.75 / 0.923 | 0.70 / 0.849 | 0.94 / 0.956 | 0.91 / 0.932 | 0.99 / 0.902 |
| f = 2.0 | 5 | 1.00 / 0.940 | 1.00 / 1.000 | 1.00 / 1.000 | 0.80 / 0.946 | 0.75 / 0.884 | 0.95 / 0.967 | 0.93 / 0.942 | 0.99 / 0.870 |

## Power, small effects (A2.2): power in the top 1% / AUROC, single pairs

| Effect | reps | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|---|
| f = 1.1 | 5 | 0.55 / 0.868 | 0.71 / 0.926 | 0.72 / 0.931 | 0.50 / 0.821 | 0.51 / 0.730 | 0.43 / 0.860 | 0.39 / 0.834 | 0.47 / 0.840 |
| f = 1.25 | 5 | 0.73 / 0.915 | 0.99 / 0.991 | 0.99 / 0.992 | 0.72 / 0.898 | 0.71 / 0.827 | 0.85 / 0.937 | 0.85 / 0.910 | 0.61 / 0.867 |
| f = 1.5 | 5 | 0.98 / 0.951 | 1.00 / 0.999 | 1.00 / 0.999 | 0.79 / 0.934 | 0.77 / 0.874 | 0.93 / 0.965 | 0.90 / 0.937 | 0.76 / 0.905 |

## Power, held-out Lung: power in the top 1% / AUROC, single pairs

| Effect | reps | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|---|
| f = 1.5 | 5 | 0.99 / 0.945 | 0.99 / 0.991 | 0.97 / 0.993 | 0.70 / 0.905 | 0.68 / 0.820 | 0.97 / 0.951 | 0.96 / 0.931 | 0.99 / 0.924 |
| f = 2.0 | 5 | 1.00 / 0.945 | 1.00 / 0.998 | 1.00 / 0.999 | 0.82 / 0.938 | 0.76 / 0.875 | 0.98 / 0.970 | 0.97 / 0.951 | 1.00 / 0.903 |

## Each tool's p inside LIANA+ on false rows (Mammary and Liver, both sides): share ≤ 0.05 / BH discoveries / false share of discoveries

| Tool | Null | runs | p ≤ 0.05 | BH discoveries | false share of discoveries |
|---|---|---|---|---|---|
| Connectome | JLP (#/B) | 12 | 0.393 [0.372, 0.415] | 72,864 | 0.940 |
| Connectome | BPT | 12 | 0.026 [0.023, 0.029] | 0 | – |
| Connectome | MDT | 12 | 0.001 [0.000, 0.002] | 2 | 0.000 |
| Connectome | LOO-t | 12 | 0.073 [0.070, 0.075] | 1,509 | 0.908 |
| NATMI | JLP (#/B) | 12 | 0.418 [0.397, 0.439] | 75,745 | 0.936 |
| NATMI | BPT | 12 | 0.029 [0.022, 0.035] | 0 | – |
| NATMI | MDT | 12 | 0.001 [0.000, 0.002] | 2 | 0.000 |
| NATMI | LOO-t | 12 | 0.072 [0.070, 0.074] | 1,498 | 0.908 |

## Real labels: calls (specificity_rank ≤ 0.05) and Jaccard with liana6

| Data set | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|
| Allen brain | 42,879 (1.00) | 39,720 (0.57) | 37,108 (0.54) | 56,037 (0.51) | 71,622 (0.41) | 46,101 (0.67) | 51,792 (0.59) | 36,815 (0.79) |
| hypothalamus | 1,456 (1.00) | 1,287 (0.69) | 1,233 (0.55) | 2,337 (0.48) | 3,033 (0.37) | 1,670 (0.52) | 2,109 (0.40) | 1,807 (0.66) |
| Liver | 33,349 (1.00) | 33,799 (0.69) | 32,741 (0.62) | 48,921 (0.56) | 60,466 (0.45) | 37,669 (0.68) | 44,598 (0.57) | 32,934 (0.81) |
| Mammary | 15,946 (1.00) | 14,284 (0.62) | 14,546 (0.55) | 21,254 (0.54) | 26,398 (0.44) | 16,605 (0.66) | 19,344 (0.58) | 12,996 (0.70) |
| TNBC | 23,466 (1.00) | 25,068 (0.61) | 27,133 (0.54) | 37,853 (0.51) | 47,111 (0.41) | 28,154 (0.64) | 33,165 (0.54) | 25,551 (0.72) |

## Biology (A2.8): spatial odds ratio, kept vs dropped [95% CI], q across all

| Variant | TNBC | Allen brain | hypothalamus |
|---|---|---|---|
| swap_null_b | 0.50 [0.29, 1.44] q 0.48 | 1.35 [0.96, 1.86] q 0.37 | 1.15 [0.72, 2.13] q 0.87 |
| swap_b2 | 1.33 [0.80, 2.22] q 0.57 | 1.48 [0.96, 2.39] q 0.32 | 1.82 [0.93, 4.33] q 0.35 |
| swap_t2 | 0.88 [0.64, 1.27] q 0.80 | 1.15 [0.86, 1.61] q 0.75 | 1.03 [0.68, 1.53] q 0.93 |
| swap3_null_a | 1.09 [0.42, 2.81] q 0.92 | 0.55 [0.37, 0.95] q 0.22 | 0.61 [0.33, 1.88] q 0.62 |
| swap3_null_b | 0.54 [0.31, 1.75] q 0.57 | 1.41 [1.06, 1.80] q 0.13 | 0.89 [0.55, 1.71] q 0.88 |
| swap3_b2 | 1.33 [0.76, 2.09] q 0.62 | 1.22 [0.95, 1.54] q 0.40 | 1.13 [0.66, 2.78] q 0.88 |
| swap3_t2 | 1.15 [0.74, 1.80] q 0.80 | 0.99 [0.84, 1.21] q 0.93 | 0.72 [0.48, 1.23] q 0.54 |
| swap4_null_a | 0.91 [0.42, 1.74] q 0.95 | 0.67 [0.50, 1.03] q 0.32 | 0.65 [0.35, 2.18] q 0.71 |
| swap4_null_b | 0.59 [0.35, 1.81] q 0.62 | 1.44 [1.05, 1.90] q 0.13 | 1.12 [0.73, 1.97] q 0.88 |
| swap4_b2 | 1.22 [0.66, 2.14] q 0.80 | 1.10 [0.90, 1.33] q 0.72 | 1.16 [0.66, 2.98] q 0.88 |
| swap4_t2 | 1.12 [0.67, 1.88] q 0.87 | 0.92 [0.78, 1.17] q 0.80 | 0.78 [0.50, 1.43] q 0.72 |

Tests with q (all) < 0.05: swap_null_b | tnbc_wu2021 progeny, swap3_null_b | tnbc_wu2021 progeny, swap3_b2 | tnbc_wu2021 progeny, swap4_b2 | tnbc_wu2021 progeny, swap_null_b | allen_brain_atlas cytosig, swap_b2 | allen_brain_atlas cytosig, swap3_null_a | allen_brain_atlas cytosig, swap3_null_b | allen_brain_atlas cytosig, swap4_null_a | allen_brain_atlas cytosig, swap4_null_b | allen_brain_atlas cytosig, swap_t2 | tabula_sapiens_mammary cytosig, swap3_null_b | tabula_sapiens_liver progeny, swap_b2 | moffitt2018_hypothalamus cytosig, swap3_null_a | moffitt2018_hypothalamus cytosig, swap3_null_b | moffitt2018_hypothalamus cytosig, swap3_b2 | moffitt2018_hypothalamus cytosig, swap4_null_a | moffitt2018_hypothalamus cytosig, swap4_b2 | moffitt2018_hypothalamus cytosig.

Downstream tests (CytoSig, PROGENy) for swap4 with bootstrap p < 0.05: swap4_null_b tnbc_wu2021 progeny (2.98, p 0.037), swap4_b2 tnbc_wu2021 progeny (5.20, p 0.005), swap4_null_a allen_brain_atlas cytosig (0.40, p 0.0005), swap4_null_b allen_brain_atlas cytosig (1.93, p 0.002), swap4_null_b tabula_sapiens_liver progeny (0.45, p 0.009), swap4_b2 tabula_sapiens_liver progeny (0.34, p 0.049), swap4_null_a moffitt2018_hypothalamus cytosig (0.04, p 0.0005), swap4_null_b moffitt2018_hypothalamus cytosig (0.08, p 0.026), swap4_b2 moffitt2018_hypothalamus cytosig (0.03, p 0.0005).

## Donors (A2.3 a): labels shuffled within donor; rows with p ≤ 0.05

Global: labels permuted over all cells (LIANA+'s scheme); donor: within each donor. MDT uses no permutations.

| Data set | Tool | reps (global / donor) | JLP: global / donor | BPT: global / donor | LOO-t: global / donor | MDT |
|---|---|---|---|---|---|---|
| Lung | CellPhoneDB | 3 / 3 | 0.236 [0.222, 0.250] / 0.064 [0.051, 0.077] | 0.060 [0.051, 0.070] / 0.004 [0.002, 0.005] | 0.062 [0.059, 0.066] / 0.062 [0.059, 0.066] | 0.006 [0.005, 0.007] |
| Lung | CellChat | 3 / 3 | 0.297 [0.273, 0.321] / 0.073 [0.053, 0.093] | 0.075 [0.060, 0.091] / 0.004 [0.003, 0.005] | 0.125 [0.121, 0.130] / 0.125 [0.121, 0.130] | 0.010 [0.009, 0.012] |
| Lung | Connectome | 3 / 3 | 0.234 [0.222, 0.246] / 0.067 [0.052, 0.082] | 0.060 [0.051, 0.070] / 0.004 [0.002, 0.005] | 0.063 [0.061, 0.065] / 0.063 [0.061, 0.065] | 0.006 [0.005, 0.007] |
| Lung | NATMI | 3 / 3 | 0.224 [0.187, 0.261] / 0.067 [0.036, 0.097] | 0.050 [0.037, 0.064] / 0.004 [-0.001, 0.008] | 0.063 [0.061, 0.065] / 0.063 [0.061, 0.065] | 0.006 [0.005, 0.007] |
| Mammary | CellPhoneDB | 3 / 3 | 0.178 [0.141, 0.216] / 0.075 [0.056, 0.094] | 0.037 [0.028, 0.046] / 0.004 [0.003, 0.005] | 0.055 [0.047, 0.064] / 0.055 [0.047, 0.064] | 0.003 [0.002, 0.005] |
| Mammary | CellChat | 3 / 3 | 0.193 [0.154, 0.232] / 0.102 [0.073, 0.130] | 0.040 [0.034, 0.046] / 0.006 [0.004, 0.007] | 0.134 [0.118, 0.150] / 0.134 [0.118, 0.150] | 0.004 [0.001, 0.007] |
| Mammary | Connectome | 3 / 3 | 0.178 [0.138, 0.218] / 0.079 [0.060, 0.099] | 0.037 [0.028, 0.046] / 0.004 [0.003, 0.005] | 0.064 [0.058, 0.069] / 0.064 [0.058, 0.069] | 0.003 [0.002, 0.005] |
| Mammary | NATMI | 3 / 3 | 0.213 [0.179, 0.248] / 0.077 [0.045, 0.108] | 0.052 [0.046, 0.058] / 0.004 [0.000, 0.008] | 0.063 [0.057, 0.070] / 0.063 [0.057, 0.070] | 0.003 [0.002, 0.005] |
| TNBC | CellPhoneDB | 3 / 3 | 0.369 [0.364, 0.374] / 0.062 [0.047, 0.076] | 0.128 [0.120, 0.136] / 0.003 [0.002, 0.005] | 0.057 [0.055, 0.060] / 0.057 [0.055, 0.060] | 0.019 [0.016, 0.023] |
| TNBC | CellChat | 3 / 3 | 0.506 [0.487, 0.524] / 0.060 [0.041, 0.078] | 0.184 [0.173, 0.196] / 0.003 [0.002, 0.004] | 0.228 [0.216, 0.240] / 0.228 [0.216, 0.240] | 0.046 [0.045, 0.047] |
| TNBC | Connectome | 3 / 3 | 0.383 [0.369, 0.396] / 0.066 [0.053, 0.080] | 0.128 [0.120, 0.136] / 0.003 [0.002, 0.005] | 0.070 [0.069, 0.071] / 0.070 [0.069, 0.071] | 0.019 [0.016, 0.023] |
| TNBC | NATMI | 3 / 3 | 0.400 [0.378, 0.422] / 0.075 [0.051, 0.098] | 0.134 [0.119, 0.149] / 0.004 [0.002, 0.007] | 0.070 [0.069, 0.071] / 0.070 [0.069, 0.071] | 0.019 [0.016, 0.023] |
| Liver | CellPhoneDB | 3 / 3 | 0.304 [0.279, 0.330] / 0.060 [0.047, 0.073] | 0.107 [0.097, 0.117] / 0.003 [0.002, 0.003] | 0.058 [0.052, 0.063] / 0.058 [0.052, 0.063] | 0.004 [0.003, 0.004] |
| Liver | CellChat | 3 / 3 | 0.377 [0.332, 0.422] / 0.075 [0.051, 0.100] | 0.146 [0.124, 0.169] / 0.004 [0.002, 0.005] | 0.123 [0.119, 0.127] / 0.123 [0.119, 0.127] | 0.005 [0.003, 0.006] |
| Liver | Connectome | 3 / 3 | 0.286 [0.258, 0.313] / 0.064 [0.053, 0.075] | 0.107 [0.097, 0.117] / 0.003 [0.002, 0.003] | 0.062 [0.060, 0.065] / 0.062 [0.060, 0.065] | 0.004 [0.003, 0.004] |
| Liver | NATMI | 3 / 3 | 0.201 [0.163, 0.239] / 0.054 [0.028, 0.079] | 0.044 [0.029, 0.060] / 0.002 [0.001, 0.003] | 0.062 [0.060, 0.065] / 0.062 [0.060, 0.065] | 0.004 [0.003, 0.004] |

## Sensitivity (A2.5)

✓: false share of the top 1% below liana6's in the same run; ✗: not below.

Calibration:

| Run | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|
| tabula_sapiens_mammary@expr_prop=0.05|receptor|0 | 0.849 | 0.631 ✓ | 0.566 ✓ | 0.761 ✓ | 0.716 ✓ | 0.846 ✓ | 0.879 ✗ | 0.860 ✗ |
| tabula_sapiens_mammary@expr_prop=0.05|ligand|0 | 0.874 | 0.719 ✓ | 0.679 ✓ | 0.852 ✓ | 0.860 ✓ | 0.939 ✗ | 0.960 ✗ | 0.939 ✗ |
| tabula_sapiens_mammary@expr_prop=0.2|receptor|0 | 0.734 | 0.489 ✓ | 0.137 ✓ | 0.496 ✓ | 0.555 ✓ | 0.830 ✗ | 0.880 ✗ | 0.849 ✗ |
| tabula_sapiens_mammary@expr_prop=0.2|ligand|0 | 0.752 | 0.387 ✓ | 0.400 ✓ | 0.638 ✓ | 0.651 ✓ | 0.925 ✗ | 0.961 ✗ | 0.917 ✗ |
| tabula_sapiens_mammary@cells=0.5|receptor|0 | 0.815 | 0.510 ✓ | 0.423 ✓ | 0.700 ✓ | 0.688 ✓ | 0.810 ✓ | 0.872 ✗ | 0.848 ✗ |
| tabula_sapiens_mammary@cells=0.25|receptor|0 | 0.771 | 0.438 ✓ | 0.371 ✓ | 0.686 ✓ | 0.633 ✓ | 0.796 ✗ | 0.859 ✗ | 0.876 ✗ |
| tabula_sapiens_mammary@labels=broad_cell_class|receptor|0 | 0.749 | 0.311 ✓ | 0.348 ✓ | 0.602 ✓ | 0.600 ✓ | 0.813 ✗ | 0.844 ✗ | 0.849 ✗ |
| tabula_sapiens_mammary@labels=broad_cell_class|ligand|0 | 0.853 | 0.580 ✓ | 0.512 ✓ | 0.809 ✓ | 0.837 ✓ | 0.923 ✗ | 0.955 ✗ | 0.957 ✗ |
| tabula_sapiens_mammary@labels=free_annotation|receptor|0 | 0.802 | 0.523 ✓ | 0.447 ✓ | 0.710 ✓ | 0.684 ✓ | 0.843 ✗ | 0.882 ✗ | 0.848 ✗ |
| tabula_sapiens_mammary@labels=free_annotation|ligand|0 | 0.839 | 0.589 ✓ | 0.525 ✓ | 0.778 ✓ | 0.792 ✓ | 0.942 ✗ | 0.963 ✗ | 0.935 ✗ |

Calibration: runs below liana6 — swap3_null_b 10 of 10, swap4_null_b 10 of 10, swap3_b2 10 of 10, swap4_b2 10 of 10, swap3_t2 2 of 10, swap4_t2 0 of 10, swap4_null_a 0 of 10.

Permutations, calibration (receptors, rep 0):

| Run | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|
| tabula_sapiens_mammary@n_perms=1000|receptor|0 | 0.804 | 0.496 ✓ | 0.379 ✓ | 0.687 ✓ | 0.656 ✓ | 0.833 ✗ | 0.871 ✗ | 0.850 ✗ |
| tabula_sapiens_mammary@n_perms=5000|receptor|0 | 0.806 | 0.503 ✓ | 0.391 ✓ | 0.687 ✓ | 0.656 ✓ | 0.833 ✗ | 0.871 ✗ | 0.841 ✗ |
| moffitt2018_hypothalamus@n_perms=1000|receptor|0 | 0.758 | 0.117 ✓ | 0.067 ✓ | 0.250 ✓ | 0.242 ✓ | 0.533 ✓ | 0.775 ✗ | 0.758 ✗ |
| moffitt2018_hypothalamus@n_perms=5000|receptor|0 | 0.758 | 0.108 ✓ | 0.067 ✓ | 0.250 ✓ | 0.242 ✓ | 0.533 ✓ | 0.775 ✗ | 0.750 ✓ |

Permutations, calibration (receptors, rep 0): runs below liana6 — swap3_null_b 4 of 4, swap4_null_b 4 of 4, swap3_b2 4 of 4, swap4_b2 4 of 4, swap3_t2 2 of 4, swap4_t2 0 of 4, swap4_null_a 1 of 4.

Power, f = 1.5 (power in the top 1% / AUROC):

| Run | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|
| tabula_sapiens_mammary@expr_prop=0.05|f2|f1.5|0 | 1.00 / 0.961 | 1.00 / 0.999 | 1.00 / 1.000 | 0.63 / 0.906 | 0.63 / 0.821 | 0.97 / 0.973 | 0.97 / 0.964 | 1.00 / 0.911 |
| tabula_sapiens_mammary@expr_prop=0.05|f2|f1.5|1 | 1.00 / 0.956 | 1.00 / 0.999 | 1.00 / 1.000 | 0.80 / 0.928 | 0.67 / 0.840 | 0.93 / 0.963 | 0.93 / 0.953 | 1.00 / 0.925 |
| tabula_sapiens_mammary@expr_prop=0.2|f2|f1.5|0 | 1.00 / 0.933 | 1.00 / 0.999 | 1.00 / 1.000 | 0.63 / 0.879 | 0.63 / 0.787 | 0.77 / 0.889 | 0.77 / 0.829 | 0.63 / 0.854 |
| tabula_sapiens_mammary@expr_prop=0.2|f2|f1.5|1 | 0.90 / 0.914 | 1.00 / 0.999 | 1.00 / 1.000 | 0.67 / 0.900 | 0.67 / 0.809 | 0.60 / 0.862 | 0.57 / 0.745 | 0.77 / 0.850 |
| tabula_sapiens_mammary@cells=0.5|f2|f1.5|0 | 0.97 / 0.949 | 1.00 / 0.996 | 1.00 / 0.997 | 0.73 / 0.920 | 0.73 / 0.853 | 0.90 / 0.932 | 0.87 / 0.908 | 0.97 / 0.911 |
| tabula_sapiens_mammary@cells=0.25|f2|f1.5|0 | 0.90 / 0.920 | 0.97 / 0.989 | 0.97 / 0.992 | 0.70 / 0.908 | 0.70 / 0.847 | 0.83 / 0.921 | 0.87 / 0.904 | 0.87 / 0.915 |

Seeds, real labels (Mammary), against seed 1337 (Spearman of specificity_rank / top-1% Jaccard):

| Run | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|
| seed1 | 0.9979 / 0.964 | 0.9936 / 0.957 | 0.9933 / 0.944 | not run | not run | 1.0000 / 1.000 | 1.0000 / 1.000 | 0.9935 / 0.928 |
| seed2 | 0.9977 / 0.979 | 0.9931 / 0.977 | 0.9934 / 0.952 | not run | not run | 1.0000 / 1.000 | 1.0000 / 1.000 | 0.9934 / 0.935 |

5,000 vs 1,000 permutations, real labels (Spearman of specificity_rank / top-1% Jaccard):

| Run | liana6 | swap3_null_b | swap4_null_b | swap3_b2 | swap4_b2 | swap3_t2 | swap4_t2 | swap4_null_a |
|---|---|---|---|---|---|---|---|---|
| real@tabula_sapiens_mammary | 0.9971 / 0.755 | 0.9953 / 0.792 | 0.9931 / 0.854 | not run | not run | 1.0000 / 1.000 | 1.0000 / 1.000 | 0.9930 / 0.805 |
| real | 0.9991 / 1.000 | 0.9961 / 0.943 | 0.9947 / 0.889 | not run | not run | 1.0000 / 1.000 | 1.0000 / 1.000 | 0.9945 / 0.962 |

"not run": A2.5's seed and permutation runs compute no MDT (`b2=False`), so the MDT variants have no second table to compare (as in AMENDMENT_7).

## Case studies (A2.9): what swap4_<n> drops among liana6's calls (s ≠ t)

Breadth: share of cell types where the gene is detected in ≥ 10% of cells; a row's breadth is the larger of its ligand's and receptor's. p: one-sided Mann–Whitney, dropped broader than kept. Overlap with swap3_<n>'s drops. Computed locally from this amendment's real tables.

| Data set | Variant | kept | dropped | breadth dropped / kept | p | dropped by both / swap4 only / swap3 only |
|---|---|---|---|---|---|---|
| TNBC | swap4_null_b | 16764 | 5518 | 0.897 / 0.586 | 5e-201 | 4869 / 649 / 1 |
| TNBC | swap4_b2 | 19455 | 2827 | 0.897 / 0.586 | 3.4e-75 | 2674 / 153 / 0 |
| TNBC | swap4_t2 | 18880 | 3402 | 0.414 / 0.724 | 1 | 3158 / 244 / 6 |
| Allen brain | swap4_null_b | 26447 | 14452 | 0.957 / 0.870 | 1.6e-121 | 12401 / 2051 / 18 |
| Allen brain | swap4_b2 | 31507 | 9392 | 0.913 / 0.913 | 1.6e-06 | 9019 / 373 / 3 |
| Allen brain | swap4_t2 | 33475 | 7424 | 0.783 / 0.913 | 1 | 6897 / 527 / 9 |
| Mammary | swap4_null_b | 9858 | 4788 | 0.875 / 0.812 | 8.1e-08 | 3999 / 789 / 51 |
| Mammary | swap4_b2 | 11896 | 2750 | 0.750 / 0.812 | 1 | 2629 / 121 / 58 |
| Mammary | swap4_t2 | 11808 | 2838 | 0.562 / 0.875 | 1 | 2714 / 124 / 146 |
| Liver | swap4_null_b | 23735 | 7638 | 0.958 / 0.917 | 7.9e-81 | 5549 / 2089 / 6 |
| Liver | swap4_b2 | 27548 | 3825 | 0.917 / 0.958 | 1 | 3533 / 292 / 0 |
| Liver | swap4_t2 | 26639 | 4734 | 0.708 / 0.958 | 1 | 4368 / 366 / 26 |
| hypothalamus | swap4_null_b | 792 | 437 | 1.000 / 0.875 | 0.12 | 291 / 146 / 0 |
| hypothalamus | swap4_b2 | 1031 | 198 | 0.812 / 1.000 | 1 | 185 / 13 / 0 |
| hypothalamus | swap4_t2 | 864 | 365 | 0.375 / 1.000 | 1 | 324 / 41 / 0 |

## Benchmark 1's tests at tool level, inside LIANA+, all four tools

p from each null, BH at 5% over the rows with a p. JLP twice: LIANA+'s p (#/B) and Benchmark 1's ((1 + #)/(B + 1)). NATMI's p are on `spec_weight`.

### Labels shuffled (any call is false): runs with any call / mean number of calls

| Data set | Tool | reps | JLP (#/B) | JLP ((1+#)/(B+1)) | BPT | MDT | LOO-t |
|---|---|---|---|---|---|---|---|
| Mammary | CellPhoneDB | 3 | 3 / 164.3 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 758.7 |
| Mammary | CellChat | 3 | 3 / 64.3 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 3,289.3 |
| Mammary | Connectome | 3 | 3 / 179.7 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 2,349.0 |
| Mammary | NATMI | 3 | 3 / 193.3 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 2,325.7 |
| Liver | CellPhoneDB | 3 | 3 / 392.7 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 2,805.3 |
| Liver | CellChat | 3 | 3 / 222.3 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 4,449.0 |
| Liver | Connectome | 3 | 3 / 458.0 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 3,606.0 |
| Liver | NATMI | 3 | 3 / 574.3 | 0 / 0.0 | 0 / 0.0 | 0 / 0.0 | 3 / 3,594.7 |

### Power (F2's design, 30 single pairs, f = 1.5, Mammary): power / FDP

| B | Tool | reps | JLP (#/B) | JLP ((1+#)/(B+1)) | BPT | MDT | LOO-t |
|---|---|---|---|---|---|---|---|
| 1,000 | CellPhoneDB | 5 | 1.00 / 0.99 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.75 / 0.18 |
| 1,000 | CellChat | 5 | 0.99 / 0.97 | 0.00 / 0.00 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.86 |
| 1,000 | Connectome | 5 | 1.00 / 0.98 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.87 |
| 1,000 | NATMI | 5 | 1.00 / 0.98 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.87 |
| 5,000 | CellPhoneDB | 5 | 1.00 / 0.98 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.75 / 0.18 |
| 5,000 | CellChat | 5 | 0.99 / 0.97 | 0.99 / 0.97 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.86 |
| 5,000 | Connectome | 5 | 1.00 / 0.98 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.87 |
| 5,000 | NATMI | 5 | 1.00 / 0.98 | 1.00 / 0.98 | 0.00 / 0.00 | 0.61 / 0.00 | 0.87 / 0.87 |

### Real labels: BH calls

| Data set | Tool | rows | JLP (#/B) | JLP ((1+#)/(B+1)) | BPT | MDT | LOO-t |
|---|---|---|---|---|---|---|---|
| Mammary | CellPhoneDB | 86,001 | 35,054 | 34,956 | 14,983 | 158 | 3 |
| Mammary | CellChat | 86,001 | 20,213 | 20,166 | 9,186 | 158 | 2,108 |
| Mammary | Connectome | 86,001 | 35,980 | 35,903 | 14,983 | 158 | 1,618 |
| Mammary | NATMI | 86,001 | 42,522 | 42,411 | 15,674 | 158 | 1,609 |
| Liver | CellPhoneDB | 196,381 | 88,904 | 88,704 | 39,270 | 9 | 296 |
| Liver | CellChat | 196,381 | 53,122 | 52,995 | 25,979 | 9 | 4,311 |
| Liver | Connectome | 196,381 | 90,109 | 89,897 | 39,270 | 9 | 3,255 |
| Liver | NATMI | 196,381 | 97,467 | 97,236 | 38,582 | 9 | 3,243 |
| TNBC | CellPhoneDB | 154,213 | 106,730 | 106,646 | 61,810 | 357 | 11 |
| TNBC | CellChat | 154,213 | 40,856 | 40,817 | 25,448 | 357 | 3,296 |
| TNBC | Connectome | 154,213 | 113,459 | 113,375 | 61,810 | 357 | 1,886 |
| TNBC | NATMI | 154,213 | 104,246 | 104,106 | 52,006 | 357 | 1,890 |
| hypothalamus | CellPhoneDB | 10,197 | 6,565 | 6,562 | 4,232 | 10 | 5 |
| hypothalamus | CellChat | 10,197 | 2,645 | 2,643 | 1,838 | 10 | 146 |
| hypothalamus | Connectome | 10,197 | 6,925 | 6,923 | 4,232 | 10 | 104 |
| hypothalamus | NATMI | 10,197 | 6,875 | 6,873 | 3,792 | 10 | 104 |
| Allen brain | CellPhoneDB | 235,366 | 106,793 | 106,524 | 35,654 | 12 | 906 |
| Allen brain | CellChat | 235,366 | 69,719 | 69,570 | 30,557 | 12 | 4,916 |
| Allen brain | Connectome | 235,366 | 100,902 | 100,645 | 35,654 | 12 | 4,538 |
| Allen brain | NATMI | 235,366 | 91,164 | 90,823 | 36,497 | 12 | 4,531 |

### Spatial truth among JLP's calls ((1 + #)/(B + 1)): odds ratio kept vs dropped [95% CI]; q: BH over all tests below (each session wrote q over its own data set; recomputed here)

| Data set | Tool | BPT | MDT | LOO-t |
|---|---|---|---|---|
| TNBC | CellPhoneDB | 0.80 [0.60, 1.11] q 0.38 | 0.41 [0.09, 1.12] q 0.29 | 4.30 [1.93, 18.07] q 0.02 |
| TNBC | CellChat | 0.76 [0.59, 1.01] q 0.29 | 0.41 [0.09, 1.23] q 0.31 | 1.12 [0.33, 1.93] q 0.93 |
| TNBC | Connectome | 0.81 [0.66, 1.04] q 0.31 | 0.41 [0.09, 1.11] q 0.29 | 0.97 [0.23, 2.02] q 0.93 |
| TNBC | NATMI | 0.79 [0.61, 1.04] q 0.29 | 0.43 [0.10, 1.13] q 0.29 | 1.02 [0.24, 2.16] q 0.97 |
| hypothalamus | CellPhoneDB | 0.89 [0.65, 1.21] q 0.61 | 3.09 [0.39, 19.74] q 0.52 | 2.05 [0.41, 20.29] q 0.57 |
| hypothalamus | CellChat | 1.04 [0.75, 1.48] q 0.93 | 3.31 [0.41, 20.25] q 0.52 | 1.03 [0.57, 1.77] q 0.93 |
| hypothalamus | Connectome | 0.84 [0.62, 1.15] q 0.51 | 3.00 [0.38, 19.18] q 0.52 | 1.27 [0.54, 2.42] q 0.65 |
| hypothalamus | NATMI | 0.99 [0.75, 1.27] q 0.93 | 3.14 [0.39, 19.74] q 0.52 | 1.33 [0.54, 2.66] q 0.61 |
| Allen brain | CellPhoneDB | 1.39 [0.98, 1.90] q 0.29 | 1.92 [0.63, 7.97] q 0.53 | 2.36 [1.32, 3.72] q 0.11 |
| Allen brain | CellChat | 1.18 [0.96, 1.44] q 0.31 | 2.22 [0.72, 8.64] q 0.46 | 1.55 [0.90, 2.40] q 0.31 |
| Allen brain | Connectome | 1.35 [0.97, 1.80] q 0.29 | 1.88 [0.60, 7.86] q 0.55 | 2.05 [1.20, 2.98] q 0.11 |
| Allen brain | NATMI | 0.94 [0.76, 1.12] q 0.68 | 1.98 [0.64, 7.99] q 0.52 | 2.21 [1.22, 3.37] q 0.11 |

Tests with q < 0.05: TNBC CellPhoneDB LOO-t.

## Checks

liana6, swap_<n> and swap3_<n> against AMENDMENT_6 (same seeds); CellPhoneDB's, CellChat's and Connectome's numbers against AMENDMENT_7; NATMI's LOO-t against Connectome's (AMENDMENT_8 expected them equal; see below).

| Check | Result |
|---|---|
| calibration (vs amendment6) | 12 of 12 runs identical |
| confirm_calibration (vs amendment6) | 6 of 6 runs identical |
| mouse_calibration (vs amendment6) | 6 of 6 runs identical |
| power (vs amendment6) | 10 of 10 runs identical |
| hard_power (vs amendment6) | 15 of 15 runs identical |
| confirm_power (vs amendment6) | 10 of 10 runs identical |
| sensitivity_calibration (vs amendment7) | 10 of 10 runs identical |
| sensitivity_power (vs amendment7) | 6 of 6 runs identical |
| sensitivity_seed (vs amendment7) | 2 of 2 runs identical |
| sensitivity_permutations (vs amendment7) | 2 of 2 runs identical |
| sensitivity_permutations_calibration (vs amendment7) | 4 of 4 runs identical |
| real calls (vs amendment6) | 5 of 5 data sets identical |
| donor_null: CellPhoneDB, CellChat, Connectome (vs amendment7) | 24 of 24 runs identical |
| tool_labels: CellPhoneDB, CellChat, Connectome (vs amendment7) | 6 of 6 runs identical |
| tool_power: CellPhoneDB, CellChat, Connectome (vs amendment7) | 10 of 10 runs identical |
| tool_real: CellPhoneDB, CellChat, Connectome (vs amendment7) | 5 of 5 runs identical |
| donor_null: NATMI's LOO-t vs Connectome's | identical in 0 of 24 runs; share of rows with p ≤ 0.01 or 0.05 differs by at most 0.0448 percentage points |
| tool_real: NATMI's LOO-t calls vs Connectome's | identical in 1 of 5 data sets; at most 12 calls apart |

AMENDMENT_8 stated that NATMI's LOO-t equals Connectome's, since within an interaction the sums over types are constant. That holds for single-gene interactions only (there the two p-values agree to rounding: LIANA+'s scores are float32). For a complex, LIANA+ uses the subunit that is lowest in each (s, t), so the sums, and with them NATMI's LOO-t, change across the interaction's rows. On the 600-cell toy data set the two p-values differ by at most 6e-8 on single-gene rows and by up to 0.22 on 32 of 80 complex rows.

