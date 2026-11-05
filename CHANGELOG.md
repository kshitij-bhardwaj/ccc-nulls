# Changelog

All notable changes to this project. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
versions follow [Semantic Versioning](https://semver.org/).

## [0.2.0] - 2026-11-04

### Added
- `ccc_nulls.rank_aggregate`: LIANA+'s `rank_aggregate` with CellChat as a sixth method and the permutation
  tools' p-values (CellPhoneDB, CellChat) from Null B, the composite null, computed from LIANA+'s own
  permutations (`ccc_nulls.liana_consensus`). Modes `swap` (recommended) and `extra`.
- Donor-aware permutations: `strata=` permutes labels within each donor, for Null A and Null B alike.
- `benchmarks/liana_final`: the pre-registered benchmark (`PLAN.md`, `AMENDMENT_1.md`, `AMENDMENT_2.md`):
  fidelity, calibration, power, biology, held-out confirmation, sensitivity, donor structure, filter
  baselines, p-value calibration, a second species, case studies, run time.
- Tabula Sapiens Lung in the dataset registry.
- `lazy_permutations`: the permutation tools score rows in blocks, so LIANA+'s per-row permutation tensor
  (2 x B x rows) is never built whole; same numbers, bounded memory (Tabula Sapiens Lung needed over 30 GB).
- Running on several machines (`--only`, `merge_results`, `RUN.md`) and on Kaggle (`kaggle/`).

### Changed
- Open Problems disabled its S3 bucket: TNBC and Allen brain come from a mirror with the original md5 sums.
- `requirements/lock-liana-final.txt`: the environment the results were produced in.
- Documentation (Sphinx), a quickstart notebook, continuous integration, BSD-3-Clause licence, `CITATION.cff`.

## [0.1.0] - 2026-10-13

### Added
- Significance layer with Null A (uniform label permutation) and the gene-stratum null, B2 (decoy genes,
  CAMERA-corrected), the composite null engine, simulations, native tool adapters, dataset registry.
