# Amendment 7: Connectome through the remaining tests

Written on 2026-11-12, before any part below had been computed, while AMENDMENT_6's last run (Lung power) was
still running.

**What had been seen when this was written:** every result of `PLAN.md` and AMENDMENTS 1–5. Of AMENDMENT_6, only
which Kaggle sessions had finished; none of its result files had been downloaded or read.

**Smoke test, disclosed.** Before this file was hashed, its code ran once on cut-down data to check that every part
finishes: Tabula Sapiens Mammary cut to 8,000 cells (donor, global, replicate 0; labels shuffled, replicate 0; power,
replicate 0) and the hypothalamus cut to 8,000 cells (real labels). Its log printed numbers, and they were seen. Two
choices below were made after seeing them, both to match Benchmark 1's conventions: JLP's calls are also counted with
the p-value (1 + #)/(B + 1) (`null_a_plus1`; LIANA+'s own p is #/B and can be 0), and tool-level power is also run with
5,000 permutations (Benchmark 1 did so for BPT, whose BH floor 1/(B + 1) is out of reach at 1,000). No number of the
smoke test is used. The full-size runs below had not been computed.

**Why.** AMENDMENT_6 put Connectome's p-value through calibration, held-out Lung, power, mouse, real data and
biology, and left out three parts: donors (A2.3), sensitivity (A2.5) and case studies (A2.9). It also checked
Connectome's own p-values only on known nulls. This amendment runs those three parts with Connectome, and the tool-level
tests that Benchmark 1 (`benchmarks/liana_nulls`) ran for CellChat and CellPhoneDB.

Code: `benchmarks/liana_final/amendment7.py`, with AMENDMENT_6's switch on: Connectome captured, T2 in every run,
and the `swap3_<n>` variants. Run: `python -m benchmarks.liana_final.amendment7`.

## Benchmark 2 parts, with Connectome

| Part | As in | Runs |
|---|---|---|
| `donor` | A2.3 (a), AMENDMENT_4, AMENDMENT_5 | labels shuffled within donor (Lung, TNBC, Mammary, Liver; 3 replicates; seeds 9000 + rep); global and within-donor permutations; share of rows with p ≤ 0.01 and ≤ 0.05 for JLP, BPT, LOO-t (both schemes) and MDT (global only: it uses no permutations), for all three tools |
| `sensitivity` | A2.5 | unchanged: seeds, `expr_prop`, cells, label resolution, 1,000 vs 5,000 permutations; `swap3_<n>` beside the earlier variants; stability against AMENDMENT_6's real tables |
| `cases` | A2.9 | among liana6's calls (s ≠ t), what `swap3_<n>` drops: count, breadth of the limiting gene (Mann–Whitney, dropped > kept), and its overlap with `swap_<n>`'s drops |

## Benchmark 1 measures for Connectome (tool level)

Connectome's own rows inside LIANA+ (`run(..., tools_only=True, connectome=True)`), p-values from each null on
`expr_prod`. Calls: BH at 5% over the tool's rows that have a p-value; JLP both with LIANA+'s p (`null_a`) and with
Benchmark 1's (`null_a_plus1`). CellChat and CellPhoneDB are recorded in the
same runs, for reference.

| Part | As in | Runs |
|---|---|---|
| `tool_labels` | Benchmark 1, labels shuffled | every label permuted (seed 4000 + rep), Mammary and Liver, 3 replicates: any call, number of calls |
| `tool_power` | Benchmark 1 power | F2's design (power.background, power.design, 30 single pairs, f = 1.5), 5 replicates, Mammary, B = 1,000 and 5,000: power and FDP (`power.metrics`) |
| `tool_real` | Benchmark 1 real data and biology | the five data sets, real labels: BH calls; among JLP's calls (`null_a_plus1`), the spatial odds ratio of those each other null keeps vs drops (A2.8's `spatial_K`, `boot_or`), TNBC, Allen brain and hypothalamus |

Known nulls at tool level (FDR at BH 5%) are AMENDMENT_6's `p_calibration` for Connectome (Mammary and Liver,
both sides, 3 replicates); they are not run again.

## Reported, not decided

Every number above is reported beside the earlier tests' numbers. No new rule is decided here: AMENDMENT_6's
decision for `swap3_<n>` stands as AMENDMENT_6 states it.

## Checks

CellChat's and CellPhoneDB's numbers in `donor` must equal AMENDMENT_2's (JLP, BPT), AMENDMENT_4's (LOO-t) and
AMENDMENT_5's (MDT) for the same tags; liana6 and the earlier variants in `sensitivity` must equal AMENDMENT_4's.
Any difference is reported.
