# Amendment 9: does the top-1% cutoff decide the result?

Written on 2026-11-16, after every result of `PLAN.md` and AMENDMENTS 1–8 was known, before any part below had been
computed.

**What had been seen when this was written:** every result of this benchmark. Since AMENDMENT_1 the consensus is
judged on its ranking: the false share of the top 1% (calibration) and power in the top 1% (planted pairs), beside
AUROC. The 1% is our choice, made in AMENDMENT_1 because `specificity_rank` is a relative rank and not a significance
threshold; it is not a convention of LIANA+. No share at any other cutoff has been computed: the result files keep the
top-1% summary only, and the rankings were not saved.

**Smoke test, disclosed.** Before this file was hashed, its code ran once on cut-down data (Mammary, 4,000 cells,
cell types with 30 cells or more, 50 permutations; one calibration run, receptors shuffled, and one planted-pair run,
f = 1.25, 2,000 cells) to check that both kinds of run finish, that the 1% values of the new measures equal the old
`false_share_top` and `power_top` (they did, for all 23 variants), and that the report is written. Read were the
variant names, the result files' keys and the report's headings. The log also printed the cut-down planted-pair
run's power and AUROC for liana6, `swap_null_b`, `extra_null_b` and `liana6_de` at 1%; those are not results of this
amendment and nothing in it uses them.

**Why.** A reader can ask whether the conclusions of AMENDMENTS 1, 6 and 8 hold only at 1%. As the cutoff grows,
every variant's false share tends to the share of false rows among all rows (the same for every variant), so gaps
must narrow; the question is whether the order changes.

Code: `benchmarks/liana_final/amendment9.py`. Run: `python -m benchmarks.liana_final.amendment9`.

## What is run

AMENDMENT_8's runs, unchanged (same code path: Connectome and NATMI captured, T2 on; same seeds, shuffles, planted
pairs and settings), so every variant of AMENDMENT_8 is present and paired:

| Part | As in AMENDMENT_8 | Runs |
|---|---|---|
| `calibration` | `calibration` | Tabula Sapiens Mammary and Liver, receptors and ligands shuffled, 3 replicates (12 runs) |
| `lung` | `lung`, calibration only | held-out Tabula Sapiens Lung, both sides, 3 replicates (6 runs) |
| `mouse` | `mouse` | Allen brain, both sides, 3 replicates (6 runs) |
| `hard_power` | `hard_power`, f = 1.25 only | A2.2's design on Mammary, 5 replicates (5 runs) |

Tool-level p-values (`p_calibration`) are not recomputed; they use no cutoff of the ranking.

## Measures added to each run

With n rows, the top fraction c is the first ⌈c·n⌉ rows of the variant's ranking (`specificity_rank` ascending,
ties by row order), exactly as `evaluate.top`.

- **Calibration**, per variant:
  - `false_share_at`: the false share of the top c, for c = 0.5%, 1%, 2%, 5% and 10%;
  - `false_share_curve`: the same at 31 fractions from 0.1% to 100%, evenly spaced on a log scale;
  - `auroc_false`: a measure with no cutoff, P(a row that is not false outranks a false row), ties one half
    (`evaluate.auroc` with the non-false rows as positives). Higher is better.
- **Power**, per variant and class of planted pair: `power_at`, the share of the planted pairs in the top c, for the
  same five c (the denominator is the planted pairs, as `power_top`); AUROC and precision at k as before.

## Reported and decided

- **Criterion 1 at each cutoff.** For each variant in `swap_<n>`, `swap3_<n>` and `swap4_<n>` (n in Null B, B2,
  T2): whether its mean false share at c is below liana6's in all four Mammary and Liver tissue × side scenarios, and
  in both Lung scenarios.
- **Verdict for BPT** (`swap_null_b`, `swap3_null_b`, `swap4_null_b`), on Mammary and Liver, then on Lung:
  - **robust to the cutoff** if criterion 1 holds at every c;
  - otherwise **holds up to c\***, the largest c up to which it holds at every smaller c, with the first c where it
    fails reported.
- **No cutoff.** Whether each variant's `auroc_false` is above liana6's in every scenario (Mammary, Liver, Lung,
  mouse); reported beside the verdict.
- Mouse and power at each c: reported, not decided.
- The report names any conclusion of the paper that changes at 2%, 5% or 10%.

## Checks

For every tag, each variant's share at 1% must equal AMENDMENT_8's `false_share_top`, and `power_at` at 1% and AUROC
must equal AMENDMENT_8's `power_top` and `auroc`. Any difference is reported first, before the results.
