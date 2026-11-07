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

(Results are filled in from `benchmarks/results/liana_final/` when the Kaggle runs finish.)
