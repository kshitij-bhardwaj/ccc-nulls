# Amendment 2: robustness, held-out confirmation and reproducibility

Written on 2026-11-04, **after** every result of `PLAN.md` and `AMENDMENT_1.md` had been seen
(`benchmarks/results/liana_final/README.md`). Seen and disclosed:
- F0 (passed on all five data sets);
- F1 and F2 as amended;
- F3;
- the verdict for both Null B variants: valid, not biologically distinguishable.

These analyses do not change those verdicts. They test how far the verdicts hold: on a tissue nobody has
looked at, on harder settings, against simpler alternatives and under other choices of parameter.
Code: `benchmarks/liana_final/amendment2.py`. Run: `python -m benchmarks.liana_final.amendment2`.

Throughout:
- **Variants and baseline** are as in `PLAN.md`; the baseline is **liana6**.
- **Measures** are as in `AMENDMENT_1.md`:
  - false share in the top 1%;
  - power in the top 1%;
  - AUROC planted vs one-sided;
  - precision at k.
- **Settings** unless stated: 1,000 permutations, seed 1337, `expr_prop` 0.1, RRA aggregation.
- **Intervals** are 95% t-intervals over replicates.

## A2.1 Held-out confirmation: Tabula Sapiens Lung

**Data and what has been seen.** Tabula Sapiens Lung (CELLxGENE, 65,847 cells, 4 donors).
- Downloaded on 2026-11-04 and not opened before this file was hashed.
- Only its catalogue entry was read: cell count, donor IDs, file size.

**The rule is frozen.** AMENDMENT_1.md's criteria 1 and 2 apply unchanged:
- the same code, seeds, scenarios and replicate counts;
- calibration: receptor-only and ligand-only shuffles, 3 replicates each;
- power: planted pairs at f = 1.5 and 2.0, 5 replicates, background and pairs drawn from Lung.

**Outcome.** For each Null B variant:
- **confirmed** if both criteria hold on Lung;
- **not confirmed** otherwise, with the failing criterion named.

**If memory runs out.** Should a full-Lung run exceed the machine's memory (8 GB), Lung is
subsampled to 30,000 cells, stratified by cell type (seed 0), and that is reported.

## A2.2 Harder power: small effects, complexes, unequal sides

**Design.**
- Same background as F2: Tabula Sapiens Mammary, 10,000 cells, labels shuffled.
- f ∈ {1.1, 1.25, 1.5}, 5 replicates.
- Each replicate plants three classes of pair, all with disjoint genes:

| Class | Pairs | What is planted |
|---|---|---|
| single | 30 | Single-gene pairs; both sides raised by f |
| complex | 10 | Pairs with a complex (≤ 3 subunits) on at least one side; every subunit raised by f |
| unequal | 10 | Single-gene pairs; ligand raised by f, receptor by 1 + (f − 1)/2 |

- Every gene is detected in ≥ 2% of cells.

**Reported.**
- Per class: power in the top 1% and AUROC planted vs one-sided.
- Overall: precision at k.
- The power gap between each variant and liana6, with its interval.

There is no pass criterion. Null B's expected power cost at small effects is what this measures.

## A2.3 Donor-aware permutations

**What changes.** LIANA+ permutes cell labels across all cells. The donor-aware version permutes
them **within each donor** (a stratified `_get_means_perms`). Both Null A and Null B are then
computed from those permutations: `null_a_donor`, `null_b_donor`.

**Donors.**
- TNBC: 10 patients (`orig.ident`).
- Tabula Sapiens Mammary: 2 donors; Liver: 4 donors; Lung: 4 donors (`donor_id`).

**(a) Known null: labels shuffled within donor.**
- Each cell type keeps its donor composition, but no type-specific expression is left within a donor.
- 3 replicates per data set.
- Measured on the tool-level p-values of CellPhoneDB and CellChat, over rows with s ≠ t and a
  positive score:
  - the share of rows with p < 0.05 and p < 0.01;
  - for Null A and Null B, each under global and under within-donor permutations.
- Expected: within-donor permutations sit closer to nominal.

**(b) Real labels.** TNBC, Mammary and Liver.
- The consensus with Null B from within-donor permutations (`swap_null_b_donor`): its top-1% overlap
  with `swap_null_b`.
- On TNBC, the F3 tests for `swap_null_b_donor`.

## A2.4 Is Null B more than a filter?

Three rankings derived from liana6:
- **`liana6_prop20` and `liana6_prop30`:** rows whose ligand (in the source) or receptor (in the
  target) is detected in fewer than 20% (30%) of cells are moved below every other row. For a
  complex, the lowest subunit counts.
- **`liana6_de`:** rows are moved below the rest unless every ligand subunit is up in the source and
  every receptor subunit is up in the target:
  - one-vs-rest t-test (scanpy `t-test_overestim_var`);
  - log fold change > 0;
  - BH-adjusted p < 0.05.

These sit beside the variants in every F1 and F2 run of this amendment (A2.1, A2.2, A2.6, A2.7).
Each Null B variant is compared head to head with each filter on false share in the top 1% and AUROC.

## A2.5 Sensitivity (Tabula Sapiens Mammary)

| Setting | Values | Runs |
|---|---|---|
| Permutations | 1,000 vs 5,000 | real labels; calibration, receptors, rep 0 |
| Seed | 1337, 1, 2 | real labels: Spearman of `specificity_rank` and Jaccard of the top 1% between seeds |
| `expr_prop` | 0.05, 0.2 (vs 0.1) | calibration, both sides, rep 0; power, f = 1.5, reps 0–1 |
| Aggregation | LIANA+'s `mean` vs `rra` | every run of this amendment, after the fact |
| Cell-type resolution | coarse `broad_cell_class`, fine `free_annotation` (vs `cell_type`) | calibration, both sides, rep 0 |
| Cells | 50% and 25% per type (seed 0) | calibration, receptors, rep 0; power, f = 1.5, rep 0 (background of 5,000 and 2,500 cells) |

Reported per setting: whether each Null B variant's false share in the top 1% stays below liana6's,
and its power and AUROC.

## A2.6 The p-values themselves

**What is rerun.** F1 on Mammary and Liver, with the same seeds, keeping every tool-level p-value.
This also checks reproducibility: the measures already reported must come out the same.

**On the false rows** (shuffled side), per tool:
- the share of p below 0.01, 0.05 and 0.1, for Null A, Null B and B2;
- the quantiles for a Q-Q plot;
- BH at 5% within each tool: the share of discoveries that are false rows, beside the share of false
  rows overall.

Null B's null holds on these rows, so its p-values should be at or above uniform. Null A's null does
not hold there.

## A2.7 A second species: Allen brain (mouse)

- F1 as amended, on the Allen brain atlas.
- Side-only genes come from LIANA+'s `mouseconsensus`.
- Receptor and ligand shuffles, 3 replicates each.
- Reported: AMENDMENT_1's criterion 1 on this data set.

## A2.8 Biology: every test, multiplicity and detectable effects

**Tests.** Every F3 test (3 spatial + 10 downstream) for every variant, plus A2.3's TNBC tests.

**p-values.** Two-sided, from the cluster bootstrap already used for the intervals (2,000
resamples):
- p = 2 · min(P(log OR* ≤ 0), P(log OR* ≥ 0)), floored at 1/2,000.

**Multiple testing.** BH within each variant's 13 tests, and across all variants and tests. No test
is primary.

**Detectable effect.** The smallest odds ratio each test could detect with 80% power at α = 0.05
(two-sided) is exp(±2.80 · SE), where SE is the bootstrap standard deviation of log OR. A test whose
detectable odds ratio is far from 1 cannot support "no difference".

## A2.9 What Null B drops, and why

For each data set, among liana6's calls with s ≠ t, compare rows swap_null_b drops with rows it keeps:
- **Limiting side:** the side with the larger Null B p (`p_L_null_b` vs `p_R_null_b`), per tool.
- **Breadth:** the share of cell types in which the limiting side's gene is detected in ≥ 10% of
  cells (lowest subunit for a complex).
- **Hypothesis:** the limiting side of dropped rows is broader than that of kept rows (one-sided
  Mann–Whitney).
- **Tables:**
  - the 20 dropped rows liana6 ranks highest;
  - the 20 rows swap_null_b calls that liana6 does not.

  Both tables list each row's side p-values and breadth. Literature claims are not generated by the
  code.

## A2.10 Reproducibility

- **Environment:** a pinned lock file of the environment the results were produced in
  (`requirements/lock-liana-final.txt`).
- **Run time:** wall time of `li.mt.rank_aggregate` against this pipeline (six methods, Null B and B2)
  on Mammary and the hypothalamus, on the same machine.
- **One command:** reproduces every result
  (`python -m benchmarks.liana_final --part all && python -m benchmarks.liana_final.amendment2 --part all`).
