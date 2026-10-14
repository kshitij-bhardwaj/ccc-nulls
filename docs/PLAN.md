# Plan

## Goal

Test whether the null model behind ligand–receptor significance changes what tools call
"communication", using each tool's own implementation and its own statistic.

## Tools (one workspace each, `tools/`)

CellPhoneDB, CellChat, NATMI, Connectome. Each is built twice:

1. **Native**: the original implementation. Its input files are produced from the dataset
   and its output is converted to the common schema (`tools/<tool>/adapter.py`).
2. **LIANA reference**: LIANA's re-implementation (`liana_reference()`, working now). It is
   the cross-check for the native version and the baseline every null is applied to.

## Steps

| # | Step | Where | Data | Status |
|---|---|---|---|---|
| 1 | Each tool's own Null A (`null_a`), on the tool's own score, with its adapter | `tools/<tool>/nulls.py`, `adapter.py` | TNBC | to do (stubs and tests in place) |
| 2 | Null A check: our Null A vs the tool's own permutation p | CellPhoneDB, CellChat | TNBC | to do; `run_tool.py` reports it |
| 3 | Each tool's own Null B (`null_b`) | `tools/<tool>/nulls.py` | TNBC | to do |
| 4 | Native adapters: dataset → native inputs → native run → common schema | `tools/<tool>/adapter.py` | TNBC | LIANA path done; native to do |
| 5 | Evaluation against ground truth (AUPRC, odds ratio, Open Problems protocol) | `benchmarks/run_tool.py` | TNBC (cytokine), mouse brain (spatial) | shared scoring in place |
| 6 | Calibration: fully permuted labels + fake LR resource | shared | Tabula Sapiens (one tissue), TNBC | TNBC done for LIANA tools; Tabula Sapiens to do |

The earlier 7-tool LIANA benchmark (`benchmarks/tnbc_openproblems.py`, shared nulls N1–N5 in
`ccc_nulls.multitool`) stays as the reference the per-tool nulls can be checked against.

## Decisions and caveats from the review of this plan

- **Two ground-truth datasets, not three.** TNBC *is* the cytokine ground truth (Open Problems
  ligand–target). The spatial truth is the mouse visual-cortex atlas (Open Problems
  source–target). Both are in `ccc_nulls.datasets`.
- **Null A exists natively in only two of the four tools.** CellPhoneDB and CellChat permute
  cluster labels. NATMI and Connectome have no permutation null on their edge score, so for
  them Null A is our implementation only, and step 4 applies to CellPhoneDB and CellChat.
- **Full permutation alone cannot calibrate.** On fully shuffled labels every null we tried sits
  near 5% (simulation S0 and real-data E5), so that test cannot separate good nulls from bad
  ones. Calibration on Tabula Sapiens should also use the fake-LR-resource test (random,
  abundance-matched non-LR gene pairs), which does separate them.
- **LIANA versions are not the originals.** LIANA re-implements the scoring (CellChat in LIANA
  uses the trimean + Hill function, but not CellChat's full pipeline), so native and LIANA results
  will differ. Each workspace documents the difference rather than forcing agreement.

## Open questions

- Stress-test metrics list (plan step 4 was cut off): current set is AUPRC, odds ratio
  (top 5%), fake-pair call rate, real/fake separation, cross-patient reproducibility,
  shuffled-label sanity. Add or drop.
- Which Tabula Sapiens tissue(s) to use for calibration (size vs coverage).

## Results so far (TNBC, LIANA versions of 7 tools)

`benchmarks/results/tnbc/`. In short:
- N1 (label permutation) calls 15–52% of fake gene pairs at 5% FDR.
- N2 (within-patient permutation) gives the same.
- N3 (abundance-matched decoys, 10,000 draws) calls 1,265–3,722 real interactions per tool
  and 0 fake ones, with a small ground-truth gain for 5 of 7 tools.
- N3 is not more reproducible across patient halves.
- N3 is ~13% at p < 0.05 on fully shuffled labels for mean-type scores (to fix).
