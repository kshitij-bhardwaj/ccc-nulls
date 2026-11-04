# AMENDMENT_2: robustness, held-out confirmation and reproducibility

Pre-registered in `benchmarks/liana_final/AMENDMENT_2.md` (committed before its runs). Run: `python -m benchmarks.liana_final.amendment2 --part all`. Measures as in AMENDMENT_1: false share in the top 1% (lower is better); power in the top 1% / AUROC planted vs one-sided / precision at k. Means over replicates [95% t-interval].

## Deviations from AMENDMENT_2

- **Permutations (A2.5).** 5,000 permutations on Tabula Sapiens Mammary need a 2 × 5,000 × 86,000 permutation tensor (about 7 GB), beyond the 8 GB machine. The comparison was run on the hypothalamus (10,197 rows), real labels and the receptor-shuffled calibration (rep 0).

