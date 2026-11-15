# Amendment 8: NATMI's specificity replaced too, through every test

Written on 2026-11-14, after every result of `PLAN.md` and AMENDMENTS 1–7 was known, before any part below had been
computed.

**What had been seen when this was written:** every result of this benchmark, Connectome's included (AMENDMENTS 6
and 7). No NATMI p-value had been computed on any data set: neither here nor under `benchmarks/liana_nulls`'s
AMENDMENT_5 (native NATMI, written by a teammate; its code is on `main` and it has not been run). The unit tests on
the 600-cell toy data set (`tests/test_natmi.py`) were run; they check identities, not results.

**Smoke test, disclosed.** Before this file was hashed, its code ran once on cut-down data to check that every part
finishes. Only that it finished, the variant names and the result files' keys were read; no number of it is used.

**Why.** LIANA+'s six methods include four with a cell-type specificity that a null can test: CellPhoneDB and
CellChat (permutation p-values, replaced since `PLAN.md`), Connectome (replaced too since AMENDMENT_6) and NATMI,
whose specificity (`spec_weight`) is not a p-value. This amendment replaces NATMI's as well and puts it through every
test Connectome went through.

Code: `ccc_nulls.liana_consensus` (`run(..., natmi=True)`, `capture_natmi`, variants `swap4_<n>`) and
`benchmarks/liana_final/amendment8.py`. Run: `python -m benchmarks.liana_final.amendment8`.

## What changes

- **NATMI's score.** LIANA+'s NATMI ranks by `expr_prod` (L_s · R_t, the same as Connectome's) and `spec_weight` =
  (L_s / Σ_k L_k) · (R_t / Σ_k R_k), the type means over all cells, Σ_k over every cell type. The nulls test
  `spec_weight`, NATMI's own score; a test of `expr_prod` would repeat Connectome's p-values. Under a permutation of
  the labels both a type's mean and the sum over types change.
- **NATMI's p-values,** from LIANA+'s permutations (the same as CellPhoneDB's: same seed, same means):
  - Null A (JLP): all labels permuted together, LIANA+'s formula #{s_b ≥ s} / B;
  - Null B (BPT): p_L with the sender side's share L / Σ L permuted and the receiver's as observed, p_R the other way
    round, p = max(p_L, p_R), (1 + k) / (B + 1); p = 1 where `spec_weight` is 0;
  - B2 (MDT): unchanged; it does not read the score, so on NATMI's rows its p equals the other tools' on the same rows;
  - T2 (LOO-t): on `spec_weight`. Within an interaction Σ_k L_k and Σ_k R_k are constant, so `spec_weight` is
    `expr_prod` times a constant and T2's p equals Connectome's. This is stated here, before any run.
- **One permutation cube.** LIANA+ builds the same mean permutations for each permuting method. With Connectome or
  NATMI captured, they are built once and shared (`shared_permutations`); CellChat's (trimean, rescaled X) are built
  as before. Same seed and arithmetic, so the same numbers (`tests/test_natmi.py`).
- LIANA+'s outputs (NATMI's `expr_prod` and `spec_weight` included) stay LIANA+'s own. So liana5, liana6 and every
  earlier variant (`swap_<n>`, `extra_<n>`, `swap3_<n>`) are unchanged.
- **New variants `swap4_<n>`,** n in Null A, Null B, B2, T2: as `swap3_<n>`, and NATMI's specificity replaced by null
  n's p (ascending). Four of LIANA+'s six methods then use null n; log2FC and SingleCellSignalR use none.
  `swap4_null_a` keeps CellPhoneDB's and CellChat's own p and gives Connectome and NATMI a Null A p.

## Parts

AMENDMENT_6's and AMENDMENT_7's code, seeds, shuffles and planted pairs, with Connectome and NATMI captured and T2
on, so `swap_<n>`, `swap3_<n>` and `swap4_<n>` come from the same runs and are paired:

| Part | As in | Runs |
|---|---|---|
| `calibration` | AMENDMENT_6 | Tabula Sapiens Mammary and Liver, receptors and ligands shuffled, 3 replicates; tool-level p of all four tools |
| `lung` | AMENDMENT_6 (A2.1) | held-out Tabula Sapiens Lung: calibration, both sides, 3 replicates; power, f = 1.5 and 2.0, 5 replicates |
| `power`, `hard_power` | AMENDMENT_6 | F2's design (f = 1.5, 2.0) and A2.2's (f = 1.1, 1.25, 1.5), 5 replicates |
| `mouse` | AMENDMENT_6 (A2.7) | Allen brain, both sides, 3 replicates |
| `real` | AMENDMENT_6 | the five data sets, real labels; tables saved |
| `biology` | AMENDMENT_6 (A2.8) | every F3 test for `swap4_<n>` beside `swap_<n>` and `swap3_<n>`; bootstrap p, BH within variant and across all |
| `donor` | AMENDMENT_7 | labels shuffled within donor (Lung, TNBC, Mammary, Liver; 3 replicates; seeds 9000 + rep); global and within-donor permutations; rows with p ≤ 0.01 and ≤ 0.05, all four tools |
| `sensitivity` | AMENDMENT_7 (A2.5) | seeds, `expr_prop`, cells, label resolution, 1,000 vs 5,000 permutations; stability against this amendment's real tables |
| `cases` | AMENDMENT_7 (A2.9) | among liana6's calls (s ≠ t), what `swap4_<n>` drops: count, breadth of the limiting gene (Mann–Whitney, dropped > kept), overlap with `swap3_<n>`'s drops |
| `tool_labels`, `tool_power`, `tool_real` | AMENDMENT_7 | Benchmark 1's tests at tool level inside LIANA+, all four tools: labels shuffled (Mammary, Liver, 3 replicates); power (F2, f = 1.5, Mammary, 5 replicates, B = 1,000 and 5,000); real labels (BH calls; spatial odds ratio among JLP's calls) |

## Reported and decided

- Per `swap4_<n>`: AMENDMENT_1's rule (criterion 1: false share of the top 1% below liana6's; criterion 2: power not
  worse), on Mammary and Liver, then on held-out Lung: **confirmed** or **not confirmed**, as for `swap3_<n>`.
- Beside each `swap4_<n>`, `swap3_<n>` and `swap_<n>` from the same runs (reported, not judged).
- Every other part: reported beside the earlier tests' numbers, not decided.

## Checks

liana6, every `swap_<n>` and every `swap3_<n>` must equal AMENDMENT_6's (calibration, Lung, power, mouse, real) and
AMENDMENT_7's (donor, sensitivity) for the same tags, and CellPhoneDB's, CellChat's and Connectome's tool-level numbers
must equal AMENDMENT_7's. Any difference is reported.

## Not run

Native NATMI (`benchmarks/liana_nulls` AMENDMENT_5: NATMI's own `ExtractEdges.py`, its fidelity check and its
Benchmark 1 stages) is a separate amendment and is not run here.
