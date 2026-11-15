# Amendment 4: T2 on the rest of AMENDMENT_2

Written on 2026-11-06, after `AMENDMENT_3.md` was hashed and before T2 had been run on any data set of this
benchmark. AMENDMENT_3 is unchanged; this file adds the parts it left out.

**What had been seen when this was written:**
- everything listed in `AMENDMENT_3.md`;
- the progress of AMENDMENT_2's Kaggle sessions, some of which (Lung calibration among them) had finished.

- a smoke test of this file's code, run before it was hashed, on 3,000 TNBC cells (stratified subsample),
  100 permutations, one within-donor shuffle (replicate 0). There, the share of rows with T2 p ≤ 0.05 was
  0.063 for CellPhoneDB and 0.263 for CellChat; Null B's was 0.045 and 0.084. The smoke test also timed T2 on the
  hypothalamus (0.2 s). Its other outputs compared against stand-in tables and carry no information about T2.

Nothing in the plan below was changed after the smoke test. No other T2 p-value had been computed on any data set.

Code: `benchmarks/liana_final/amendment4.py`. Run: `python -m benchmarks.liana_final.amendment4`.

## The null and how it is applied

As in `AMENDMENT_3.md`:
- T2 from PR #11, unchanged;
- applied to CellPhoneDB's `lr_means` and CellChat's `lr_probs` over LIANA+'s scored rows;
- `min_ref` = 3, untestable rows get p = 1;
- variants `swap_t2` and `extra_t2`.

**T2 uses no permutations.** So `swap_t2`, whose two permutation tools are replaced by T2, does not depend on
the permutation seed, the number of permutations or the permutation scheme. `extra_t2` keeps the tools' own
p-values and does. Where a part varies only the permutations, `swap_t2` is still reported, as a check that it
is unchanged.

Every run reuses AMENDMENT_2's code, seeds, shuffles and planted pairs, with T2 switched on. So liana6 and the
Null B variants come out of the same runs. Settings: 1,000 permutations, seed 1337, `expr_prop` 0.1, RRA.

## Parts

| Part | Mirrors | Runs | Reported |
|---|---|---|---|
| `lung` | A2.1 | Tabula Sapiens Lung, full: calibration, both sides, 3 replicates; power, f = 1.5 and 2.0, 5 replicates | **confirmed** if AMENDMENT_1's criteria 1 and 2 hold on Lung, else **not confirmed** with the failing criterion named; per T2 variant |
| `donor` (a) | A2.3 (a) | labels shuffled within donor, TNBC, Mammary, Liver and Lung, 3 replicates, global permutations | per tool, rows with s ≠ t and a positive score: share with p < 0.05 and p < 0.01 for T2, beside Null A and Null B |
| `donor` (b) | A2.3 (b) | real labels, TNBC, Mammary and Liver, within-donor permutations | top-1% and call Jaccard of `swap_t2` and `extra_t2` against the same variant under global permutations (AMENDMENT_3's `real` tables) |
| `sensitivity` | A2.5 | Mammary: seeds 1 and 2; `expr_prop` 0.05 and 0.2; 50% and 25% of cells; coarse and fine labels; 1,000 vs 5,000 permutations (hypothalamus, and Mammary) | per setting, whether each T2 variant's false share in the top 1% stays below liana6's, and its power and AUROC; for seeds and permutations, Spearman of `specificity_rank` and top-1% Jaccard against the seed-1337, 1,000-permutation run |
| `mouse` | A2.7 | Allen brain, receptor and ligand shuffles, 3 replicates | AMENDMENT_1's criterion 1, per T2 variant |
| `cases` | A2.9 | the five data sets, AMENDMENT_3's `real` tables | see below |
| `timing` | A2.10 | Mammary and the hypothalamus | wall time of T2 on the captured tools (median of 3), beside A2.10's times |

**`cases` for T2.** T2 tests the row as a whole, not one side, so A2.9's limiting side does not apply. Among
liana6's calls with s ≠ t, rows `swap_t2` drops are compared with rows it keeps:
- **Breadth:** for the ligand and for the receptor, the share of cell types in which the gene is detected in
  ≥ 10% of cells (lowest subunit for a complex); the row's breadth is the larger of the two.
- **Hypothesis:** dropped rows are broader than kept rows (one-sided Mann–Whitney).
- **Overlap with Null B:** rows dropped by both `swap_t2` and `swap_null_b`, by T2 only and by Null B only,
  and the Jaccard of the two dropped sets.
- **Tables:** the 20 dropped rows liana6 ranks highest, and the 20 rows `swap_t2` calls that liana6 does not,
  with each row's T2 p-value per tool and its breadths.

## Checks

The liana6 and Null B columns repeat AMENDMENT_2's runs and must equal them:
- `lung`: `amendment2/confirm_calibration.json` and `confirm_power.json`;
- `donor` (a): Null A and Null B under global permutations, `amendment2/donor_null.json`;
- `mouse`: `amendment2/mouse_calibration.json`;
- `sensitivity`: `amendment2/sensitivity_*.json`.

Any difference is reported.

## Not run

- the F3 tests of the donor-aware variants on TNBC (A2.3 (b), last item);
- A2.8's multiplicity and detectable-effect analysis for the T2 variants.
