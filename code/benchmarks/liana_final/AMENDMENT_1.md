# Amendment 1 to the final pipeline's plan

Written after some results were already seen. They are disclosed here:

- **F0 and the calls on all five data sets.** F0 passed everywhere.
- **F3 biology, complete.** The point estimates were in the log. This amendment leaves F3 unchanged.
- **One F2 power replicate (rep 0, f = 1.5 and 2.0).** Two flaws in the plan showed up:
  1. **Only 20 of the 30 planted pairs are rows of LIANA+'s consensus.** They were picked from
     CellChatDB ∩ CellPhoneDB, as in `benchmarks/liana_nulls`. Power was therefore capped at 0.67,
     and criterion 2 (power ≥ 0.8) could not be met by construction.
  2. **`specificity_rank ≤ 0.05` is not a significance threshold.** The RRA rank is relative. On data
     with shuffled labels, where only the 33 planted rows are real, every variant called about 17% of
     all rows (≈ 17,400 of 103,439). FDR was ≈ 1.0 for all of them, so the rule cannot tell them apart.

The F1 calibration stage was stopped before it ran. The project lead chose the changes below.

## Changes

**F2 power.**
- The 30 pairs are drawn from LIANA+'s consensus resource: single-gene pairs whose two genes are each
  detected in ≥ 2% of cells. Otherwise the design is the same: disjoint genes, random source and target,
  the same planting and the same seeds.
- The rows present are reported.
- The variants are judged on their **ranking** (`specificity_rank`, ascending), with no threshold:
  - **power in the top 1%**: share of the planted rows ranked in the top 1% of all rows;
  - **AUROC (planted vs one-sided)**: probability that a planted row outranks a one-sided row (one side
    raised, the other not). This is the specific claim Null B makes.
  - **precision at k**: share of true rows among the top k, where k = number of true rows.

**F1 calibration.**
- The ligand-only and receptor-only scenarios are kept: 2 tissues × 2 sides × 3 replicates.
- The measure becomes the **false share in the top 1%**: share of rows in the top 1% of the ranking
  whose shuffled side has every subunit shuffled. The share among all rows is reported beside it.
- The shuffled-labels scenario is dropped. Every row is null there, so a ranking measure has no target.

**F3 biology.** Unchanged: calls at `specificity_rank ≤ 0.05`, kept vs dropped among the baseline's calls.

## Decision rule (replaces the plan's)

For `swap_null_b` and `extra_null_b` separately:

1. **Calibration:** the false share in the top 1% is below liana6's in all four tissue × side scenarios
   (mean over replicates).
2. **Power:** at f = 2.0, power in the top 1% ≥ 0.8, and AUROC (planted vs one-sided) ≥ liana6's.
3. **Biology:** as in the plan.

The verdicts are as in the plan.
