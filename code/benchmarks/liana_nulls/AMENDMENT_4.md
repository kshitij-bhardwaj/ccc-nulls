# Amendment 4: the leave-one-out t-test through every stage

Written on 2026-11-09, after every result of `PLAN.md` and Amendments 1–3 was known, and before the leave-one-out
t-test had been computed on any data set of this benchmark.

**What had been seen when this was written:**
- every result of this benchmark (`benchmarks/results/liana_nulls/README.md`);
- the leave-one-out t-test's results in the other benchmark (`benchmarks/liana_final`, AMENDMENTS 3 and 4);
- the run time of the tools' scoring on Tabula Sapiens Mammary (12 s CellPhoneDB, 38 s CellChat), measured to plan
  this amendment; no p-value was computed.

**Why.** Every other null of this benchmark went through calibration, power, real data and biology. The
leave-one-out t-test (code name T2, from PR #11; called LOO-t here, to avoid confusion with this plan's power stage,
also called T2) did not. This amendment runs it through the same stages. It also adds the ligand-shuffled
calibration for Null B2, which Amendment 3 ran with receptors only.

Code: `benchmarks/liana_nulls/amendment4.py`. Run: `python -m benchmarks.liana_nulls.amendment4`.

## The test

- `ccc_nulls.nullb_type_restd.p_values` on the tool's score (`score`), unchanged from PR #11: for each interaction,
  a t-statistic of the (source, target) row against the interaction's other tested rows, leaving the row out.
- Tested rows: the tool's BH family (`family`), the rows Null A and Null B are tested on. Rows with fewer than three
  reference rows get p = 1.
- Calls: BH at α = 0.05 over the family, as for Null A and Null B.
- No permutations, so the tools are run with B = 10 only to obtain their scores; the score does not depend on B.

## Stages

| Stage | As in | Runs | Reported |
|---|---|---|---|
| calibration | Amendment 3's validity harness (A3.3) | Tabula Sapiens Mammary and Liver, both tools, 10 replicates each of: receptor genes shuffled (seeds 5000 + rep, A3.3's), ligand genes shuffled (ligand-only genes, seeds 6000 + rep), labels shuffled (seeds 4000 + rep, A3.3's) | FDR at BH 0.05 on rows whose shuffled side is null (mean, 95% t-CI); with labels shuffled, P(any call). LOO-t and Null B2 on the same shuffles |
| power | the plan's T2 stage | the plan's planted data: replicates 0–4, f = 1.5, Tabula Sapiens Mammary | power and FDR, BH at 0.05 (mean, 95% t-CI) |
| real data | the plan's T0 | the five data sets, both tools, the saved core rows | calls; calls outside Null A's |
| biology | the plan's T3 primary and Amendment 1 (A1) | TNBC, Allen brain, hypothalamus; both tools | odds ratio of colocalisation, Null A's calls kept vs dropped by LOO-t, 2,000 bootstrap resamples of (s, t) |

Null A's and Null B's calibration stays as the plan computed it (20 replicates, the tools' calibration scripts); the
harness above uses 10 replicates and other seeds. Both shuffle one side's genes across cells, so the false rows are
defined the same way.

## Decision

The plan's rule, applied to LOO-t as written for Null B: **valid** if the FDR's upper 95% bound is ≤ 0.05 in every
calibration setting and in the power stage; the biology criteria as in the plan.
