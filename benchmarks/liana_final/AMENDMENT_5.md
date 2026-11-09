# Amendment 5: the same tests for every null

Written on 2026-11-09, after every result of `PLAN.md` and AMENDMENTS 1–4 was known, before any part below had been
computed.

**What had been seen when this was written:** every result of this benchmark, including T2's in AMENDMENTS 3 and 4,
and Null B2's in AMENDMENT_2. No result of the parts below had been computed.

**Why.** Null B went through every test of AMENDMENT_2. T2 (PR #11's leave-one-out t-test) and B2 (decoy genes,
CAMERA-corrected) each missed some:
- T2: A2.8's multiplicity and detectable-effect analysis, and the TNBC test of its donor-aware variant (AMENDMENT_4,
  "Not run");
- B2: the donor test A2.3 (a), and A2.9's "what it drops".

This amendment runs exactly those, so that every null has the same tests.

Code: `benchmarks/liana_final/amendment5.py`. Run: `python -m benchmarks.liana_final.amendment5`.

## Parts

| Part | Mirrors | Runs | Reported |
|---|---|---|---|
| `biology` | A2.8 | every F3 test of AMENDMENT_2, now with `swap_t2` and `extra_t2` (AMENDMENT_3's real tables) and `extra_t2_donor` on TNBC (AMENDMENT_4's within-donor table) | odds ratio, bootstrap p, BH q within each variant and across all tests, detectable odds ratios |
| `donor_b2` | A2.3 (a) | labels shuffled within donor (A2.3's shuffles, seeds 9000 + rep), TNBC, Mammary, Liver, Lung, 3 replicates | per tool, rows with s ≠ t and a positive score: share with B2 p < 0.05 and < 0.01, among all rows and among rows B2 tests |
| `cases_b2` | A2.9 / AMENDMENT_4's `cases` | the five data sets, AMENDMENT_3's real tables | rows `swap_b2` drops among liana6's calls: breadth of the larger side, one-sided Mann–Whitney (dropped broader), overlap with `swap_null_b`'s drops |

Details:
- **`biology`:** across-all BH covers AMENDMENT_2's tests plus the T2 tests, so its q-values differ from
  `amendment2/biology.json`'s; the within-variant q-values of AMENDMENT_2's variants do not change. The odds ratios
  and p-values of AMENDMENT_2's variants must equal `amendment2/biology.json`'s (same tables, same bootstrap seed);
  any difference is reported. `swap_t2` does not depend on the permutation scheme, so only `extra_t2_donor` is
  added.
- **`donor_b2`:** B2 uses no permutations, so there is one scheme. LIANA+ runs with 10 permutations, which change
  no score.

## Not run

Nothing of AMENDMENT_2's tests is left out for Null B, B2 or T2 after this amendment, except tests that do not
apply: global versus within-donor permutations for the two permutation-free nulls, and the filters (A2.4), which are
not nulls.
