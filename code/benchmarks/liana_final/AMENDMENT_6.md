# Amendment 6: Connectome's p-value replaced too

Written on 2026-11-11, after every result of `PLAN.md` and AMENDMENTS 1–5 was known, before any part below had been
computed.

**What had been seen when this was written:** every result of this benchmark, and PR #13's proposal (Connectome as a
permutation method inside the consensus). No Connectome p-value had been computed on any data set. The unit tests on
a 600-cell toy data set (`tests/test_connectome.py`) were run; they check identities, not results.
A smoke test of this file's code ran once before it was hashed: Tabula Sapiens Mammary cut to 8,000 cells, receptors
shuffled, replicate 0, written to a scratch folder. Only that it finished, the variant names and the result files'
keys were read; none of its numbers was looked at, and none is used.

**Why.** LIANA+'s six methods include three with a cell-type specificity that a null can test: CellPhoneDB and CellChat
(permutation p-values, replaced since `PLAN.md`) and Connectome, whose specificity (`scaled_weight`, a mean z-score)
is not a p-value. This amendment asks whether replacing Connectome's specificity as well helps the consensus.

Code: `ccc_nulls.liana_consensus` (`run(..., connectome=True)`, `capture_connectome`, variants `swap3_<n>`) and
`benchmarks/liana_final/amendment6.py`. Run: `python -m benchmarks.liana_final.amendment6`.

## What changes

- LIANA+ permutes for Connectome as it does for CellPhoneDB (same seed, same means), and Connectome's two scores
  (`expr_prod`, `scaled_weight`) stay LIANA+'s own. So liana5, liana6 and every earlier variant are unchanged
  (checked: `tests/test_connectome.py`, and below).
- On Connectome's score `expr_prod`, every null of this benchmark gives a p-value: Null A (all labels permuted, LIANA+'s
  formula #{s_b ≥ s} / B), Null B (p = max(p_L, p_R)), B2 (decoy genes, as for the other tools: it does not depend on
  the score) and T2 (leave-one-out t on `expr_prod`).
- New variants `swap3_<n>`, n in Null A, Null B, B2, T2: as `swap_<n>`, and Connectome's specificity replaced by
  null n's p (ascending). `swap3_null_a` keeps CellPhoneDB's and CellChat's own p and gives Connectome a Null A p.

## Parts

AMENDMENT_2's and AMENDMENT_3's code, seeds, shuffles and planted pairs, with Connectome captured and T2 on:

| Part | As in | Runs |
|---|---|---|
| `calibration` | AMENDMENT_3 (F1 + A2.6) | Tabula Sapiens Mammary and Liver, receptors and ligands shuffled, 3 replicates; tool-level p of all three tools |
| `lung` | A2.1 | held-out Tabula Sapiens Lung, full: calibration, both sides, 3 replicates; power, f = 1.5 and 2.0, 5 replicates |
| `power`, `hard_power` | AMENDMENT_3 | F2's design (f = 1.5, 2.0) and A2.2's (f = 1.1, 1.25, 1.5), 5 replicates |
| `mouse` | A2.7 | Allen brain, both sides, 3 replicates |
| `real` | AMENDMENT_3 | the five data sets, real labels; tables saved |
| `biology` | A2.8 | every F3 test for the `swap3_<n>` variants beside AMENDMENT_5's; bootstrap p, BH within variant and across all |

## Reported and decided

- Per `swap3_<n>`: AMENDMENT_1's rule (criterion 1: calibration, false share of the top 1% below liana6's; criterion 2:
  power not worse), on Mammary and Liver, then on held-out Lung: **confirmed** or **not confirmed**, as for Null B.
- Beside each `swap3_<n>`, `swap_<n>` from the same runs, and the paired difference (reported, not judged).
- Tool level: Connectome's p on false rows (share ≤ 0.05; BH discoveries), per null.

## Checks

liana6 and every earlier variant in these runs must equal AMENDMENT_3's and AMENDMENT_4's (same seeds); any
difference is reported.

## Not run

Donor structure (A2.3), sensitivity (A2.5) and case studies (A2.9) for the `swap3` variants.
