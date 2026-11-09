# Reproducing the benchmark

Every number in [Benchmark](benchmark.md) comes from a JSON file under `benchmarks/results/liana_final/`, written
by the code in `benchmarks/liana_final/`. Each file records the SHA-256 of the plan and amendments it was produced
under.

## Pre-registration

| File | Fixed | What it covers |
|---|---|---|
| `benchmarks/liana_final/PLAN.md` | before any run | design, F0 fidelity, F1 calibration, F2 power, F3 biology, decision rule |
| `benchmarks/liana_final/AMENDMENT_1.md` | after one power replicate (disclosed) | F1 and F2 judged on rankings |
| `benchmarks/liana_final/AMENDMENT_2.md` | after all of the above (disclosed), before its own runs | held-out tissue, sensitivity, donor structure, filters, p-value calibration, mouse, multiplicity, case studies, run time |

Deviations from an amendment are listed in the results README beside the result they affect.

## Commands

```bash
# environment (the one that produced the results)
pip install -r requirements/lock-liana-final.txt && pip install -e . --no-deps

# data: ~9 GB in data/ (or $CCC_NULLS_DATA), md5-checked
python -m ccc_nulls.datasets tnbc_wu2021 allen_brain_atlas tabula_sapiens_mammary tabula_sapiens_liver tabula_sapiens_lung
# the hypothalamus (scRNA GSE113576 + MERFISH): see data/README.md

# plan + AMENDMENT_1 (F0-F3), then AMENDMENT_2
python -m benchmarks.liana_final --part all
python -m benchmarks.liana_final.amendment2 --part all
```

Or `make reproduce`. On an 8 GB, 8-core laptop the plan takes about 10 hours and AMENDMENT_2 about a day; every
part resumes where it stopped.

## Compute

| | |
|---|---|
| Machine | Apple silicon laptop, 8 cores, 8 GB RAM |
| Permutations | 1,000 (LIANA+'s default), seed 1337 |
| Run time | see [Benchmark](benchmark.md) |
