# Amendment 3: a second Null B candidate, T2, in the same pipeline

Written on 2026-11-06, before T2 had been run on any data set of this benchmark.

**What had been seen when this was written:**
- every result of `PLAN.md` and `AMENDMENT_1.md`;
- the AMENDMENT_2 results computed on a laptop (`amendment2/README.md`);
- the Kaggle re-runs of the five data sets and of F1 on Mammary, which matched the laptop exactly.

**What had been read about T2:** only its code and docstrings. T2's own simulation results, kept on its author's
branch, were not read.

Code: `benchmarks/liana_final/amendment3.py`. Run: `python -m benchmarks.liana_final.amendment3`.

## The null

**T2, the restandardised cell-type-ID permutation.**
- **Source:** PR #11, commit `a0a843a`. It is taken unchanged: `src/ccc_nulls/type_perm.py` and
  `src/ccc_nulls/nullb_type_restd.py`.
- **The test:** for one interaction, the score at the pair (s, t) is compared with the same interaction's scores at
  every other cell-type pair. This is a leave-one-out t statistic, referred to t with m − 1 degrees of freedom.
- **What it uses:** it needs no permutations. It uses only the tool's score at every pair.

**Applied here, as PLAN.md applies its other nulls:** only the p-value of the two permutation tools is replaced.
- **Scores:** CellPhoneDB uses LIANA+'s `lr_means`; CellChat uses `lr_probs`.
- **Reference set:** the interaction's other (source, target) rows that LIANA+ scored, i.e. those that pass
  `expr_prop` = 0.1. This is T2's "tested" set.
- **Minimum reference size:** `min_ref` = 3, T2's default. A row with fewer than 3 references cannot be tested
  and gets p = 1.
- **Variants:** `swap_t2` and `extra_t2`, built by `liana_consensus.variants` exactly as `swap_null_b` and
  `extra_null_b` are.

## Runs

Every run reuses AMENDMENT_2's code: the same seeds, shuffles, backgrounds and planted pairs. So liana6 and the
Null B variants come out of the same runs as T2.

Settings: 1,000 permutations, seed 1337, `expr_prop` 0.1, RRA aggregation.

| Part | Runs | Measures |
|---|---|---|
| `calibration` | Tabula Sapiens Mammary and Liver, receptor-only and ligand-only genes shuffled, 3 replicates (F1's and A2.6's seeds) | false share in the top 1% (AMENDMENT_1); per tool, share of false rows with p ≤ 0.01, 0.05 and 0.1, and BH discoveries at 5% (A2.6) |
| `power` | F2's design: Mammary, shuffled labels, 30 planted pairs, f = 1.5 and 2.0, 5 replicates | power in the top 1%, AUROC planted vs one-sided, precision at k |
| `hard_power` | A2.2's design: f = 1.1, 1.25 and 1.5, 5 replicates; single, complex and unequal pairs | as `power`, per class |
| `real` | the five data sets of PLAN.md | calls per variant, Jaccard with liana6, T2's untestable share |
| `biology` | F3 on the `real` tables | spatial odds ratio of kept vs dropped liana6 calls (95% bootstrap CI), pair-level AUPRC, downstream odds ratio |

## Decision rule

AMENDMENT_1.md's rule, unchanged, is applied to `swap_t2` and `extra_t2` separately:

1. **Calibration:** the mean false share in the top 1% is below liana6's in all four tissue × side scenarios.
2. **Power:** at f = 2.0, power in the top 1% is ≥ 0.8, and AUROC planted vs one-sided is ≥ liana6's.
3. **Biology:** the spatial odds-ratio CI lies above 1 in at least one data set and below 1 in none.

The verdicts are as in PLAN.md:
- "improves LIANA+" if 1–3 all hold;
- "valid, not biologically distinguishable" if 1 and 2 hold;
- "not supported" otherwise.

**T2 against Null B is reported, not judged.** On the same runs, the report gives, per measure, the
difference T2 − Null B for swap against swap and for extra against extra. Each difference is shown with its mean
and 95% t-interval over replicates.

## Checks

The runs repeat earlier runs with the same seeds, so the columns for liana6, `swap_null_b` and `extra_null_b`
must equal the earlier results:
- **calibration:** A2.6, `amendment2/recalibration.json`;
- **hard_power:** A2.2, `amendment2/hard_power.json`;
- **real:** the tables of the five `run` parts, compared rank by rank;
- **power:** F2, `power.json`. These are compared as means, because F2 stored summaries.

Any difference is reported. The comparisons with T2 use this amendment's own runs only.

## Not run

T2 is not run on these parts of AMENDMENT_2:
- Lung (A2.1);
- donor-aware permutations (A2.3);
- sensitivity (A2.5);
- mouse (A2.7);
- the case studies (A2.9);
- timing (A2.10).

A later amendment can add them.
