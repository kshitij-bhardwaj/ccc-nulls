# Earlier runs replaced by the Kaggle runs

The plan's result files were first produced locally (commit b74aa63, 2026-11-03) and committed. They were rerun on
Kaggle from commit 17272a8 (`kaggle/README.md`), and the Kaggle files replace them, so that every result comes
from one commit in the locked environment.

Apart from wall times, the write date and the SHA-256 of AMENDMENT_1 (which the later code records), the two runs
agree everywhere except in `biology.json`, at these entries (bootstrap over (s, t) pairs, 2,000 resamples):

| Entry | Local, b74aa63 | Kaggle, 17272a8 |
|---|---|---|
| allen_brain_atlas, `swap_null_b`, odds ratio, upper 95% limit | 7.279 | 7.266 |
| moffitt2018_hypothalamus, `extra_null_a`, odds ratio, upper 95% limit | 10.62 | 9.50 |
| moffitt2018_hypothalamus, `extra_null_a`, resamples with no finite odds ratio | 6 | 7 |
| moffitt2018_hypothalamus, `extra_null_b`, odds ratio | 3.33 | 0.00 |
| moffitt2018_hypothalamus, `extra_null_b`, resamples with no finite odds ratio | 13 | 11 |
| moffitt2018_hypothalamus, `swap_b2`, resamples with no finite odds ratio | 22 | 25 |
| tnbc_wu2021, `swap_b2`, odds ratio, upper 95% limit | 5.135 | 5.147 |

**Why they differ:** between the two commits, the permutation tools were changed to score rows in blocks
(commit 36ca216; Lung needed more than 30 GB otherwise). That changes the order of floating-point sums, so a few
p-values at the 0.05 cut move across it, and with them a few calls. The hypothalamus odds ratios rest on a handful
of calls (up to a quarter of resamples had none finite), so one call moves the estimate from 3.33 to 0.
No conclusion in `docs/benchmark.md` depends on these entries.

## Runs on the IITK lab server

AMENDMENT_2's and AMENDMENT_4's Mammary runs with 5,000 permutations needed more memory than a Kaggle session
has, and ran on a lab server (2 cores; `machine`: "IITK lab server, 2 cores"). There, AMENDMENT_3's real-label
Mammary run had to be repeated first. That repeat agreed with the Kaggle run on every value except
`same_as_run`, which was false on the server. The server compared against the plan's local tables (commit b74aa63,
copied with the data), not the Kaggle tables, and those differ as described above. The Kaggle run is kept.
