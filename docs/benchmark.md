# Benchmark

Two pre-registered benchmarks, presented in this order. Each has its own `PLAN.md` and amendments; the amendment
numbers of one say nothing about the other's.

## 1. Each tool's own p-values (`benchmarks/liana_nulls`)

The question: is the p-value of CellPhoneDB and CellChat valid under Null A and under Null B?

- **Fidelity.** Three columns per tool: the tool itself (CellChat 2.1.2 in R, CellPhoneDB 5.0.1), our Python Null A
  on the tool's own permutations, and Null B. Columns 1 and 2 give the same p on every row of all five data sets.
- **Calibration.** Tabula Sapiens Mammary and Liver, one side's genes shuffled (a known null), BH at 5%:

  | Tool | Null A's FDR | Null B's FDR |
  |---|---|---|
  | CellChat | 0.27–0.75 | 0.00–0.03 |
  | CellPhoneDB | 0.96–1.00 | 0.00 |

- **Power.** With BH, Null B is limited by the number of permutations: CellChat at f = 1.5 has power 0.19 at
  B = 1,000 and 0.97 at B = 5,000 (FDR 0.02).
- **Biology.** Null B is not distinguishable from Null A on the spatial and downstream-response truths.
- **Also compared:** B2 (decoy genes, CAMERA-corrected) and PR #9's CellPhoneDB Null B (its AMENDMENT_3).
- **T2 (AMENDMENT_4).** Through every stage: FDR on known nulls 0.20–0.99, calls whenever labels are shuffled,
  power 0.51–0.58 at FDR 0.24–0.39. Not valid. B2 with ligands shuffled: FDR 0.

Verdict under its plan: **valid, not biologically distinguishable from Null A.** Full results:
`benchmarks/results/liana_nulls/README.md`.

## 2. LIANA+'s consensus ranking (`benchmarks/liana_final`)

The question: with Null B inside LIANA+'s consensus, does the ranking improve?

Variants:
- **liana6:** LIANA+'s five methods plus CellChat;
- **swap:** CellPhoneDB's and CellChat's p replaced by the null's;
- **extra:** the null's p added as a 7th and 8th method.

Nulls: Null B, B2 and T2 (PR #11; AMENDMENTS 3 and 4).

Parts:
- **PLAN.md and AMENDMENT_1:** calibration, power and biology.
- **AMENDMENT_2:**
  - held-out Lung;
  - harder power;
  - donor-aware permutations;
  - simple filters;
  - sensitivity;
  - tool-level p-values;
  - mouse.

- **AMENDMENTS 3 and 4:** T2 on every part above.

Results (false share of the top 1% on known nulls, lower is better):

| | liana6 | swap: Null B | 7th/8th: Null B | swap: T2 |
|---|---|---|---|---|
| Tabula Sapiens Mammary and Liver | 0.80–0.86 | 0.60–0.69 | 0.69–0.73 | 0.80–0.89 |
| Held-out Tabula Sapiens Lung | 0.83–0.88 | 0.68 | 0.73–0.75 | 0.87–0.89 |
| Allen brain (mouse) | 0.84–0.93 | 0.73–0.83 | 0.69–0.81 | 0.87–0.94 |

- **Confirmation.** On held-out Lung both Null B variants meet AMENDMENT_1's criteria (calibration and power);
  the T2 variants fail calibration.
- **Power.** Null B ranks planted pairs higher: power in the top 1% at f = 1.25 is 0.99 against 0.73, and on Lung at
  f = 2 the AUROC is 0.997 against 0.945.
- **Tool-level p inside LIANA+.** On false rows, LIANA+'s own p is ≤ 0.05 for 38% (CellPhoneDB) and 19% (CellChat);
  Null B's for 2.6% and 1.3%. With BH at 1,000 permutations Null B makes no discoveries; it needs about 5,000.
- **Filters.** A detection filter changes nothing. A both-sides DE filter removes more false rows (18 of 18
  settings) but ranks planted pairs lower than Null B (AUROC, 9 of 9).
- **Donors.** With labels shuffled within donor, global permutations give LIANA+'s p ≤ 0.05 on 18–37% of rows
  (CellPhoneDB) and 19–51% (CellChat); within-donor permutations give 6–10% for Null A and ≤ 0.6% for Null B.
  `rank_aggregate` now permutes within donor by default.
- **Sensitivity.** Null B beats the baseline in all 16 runs (seeds, `expr_prop`, subsampling, labels, 5,000
  permutations on Mammary and the hypothalamus, mouse).
- **Biology.** No spatial test is significant after BH for any variant (19 tests); detectable odds ratios are about
  1.6–3. Null B drops rows whose limiting gene is broadly expressed.
- **T2.** About a thousand times faster (2 s per data set), permutation-free, and near nominal for CellPhoneDB under
  donor structure (5.5–6.2%), but not calibrated in the consensus: it beats the baseline in 4 of 14 sensitivity runs.
- **Ceiling.** Four of the six methods use no null, so 60–83% of the top 1% stays false even with Null B.

- **AMENDMENT_5.** B2 and T2 get the tests only Null B had: T2's biology with multiplicity correction (no spatial
  test significant), B2 under within-donor shuffles (0.3–4.6% at p ≤ 0.05), and what B2 drops (broadly
  expressed genes in 4 of 5 data sets).

Verdict under the plan: **Null B is valid and gives a better ranking in every setting tested; it is not
biologically distinguishable on the available truths. T2: not supported.** Full results:
`benchmarks/results/liana_final/README.md` and `amendment2/` to `amendment4/`; reruns: `RERUNS.md`.
