# Null B on CellPhoneDB (shared composite null through the CellPhoneDB adapter)

Script: `benchmarks/null_b_cellphonedb.py`. CellPhoneDB 5.0.1, database v5.0.0, B = 1,000, BH at
α = 0.05 over the adapter's tested rows (CellPhoneDB's expression gate) for both nulls.

- Null A is the CellPhoneDB adapter's own (`cpdb_null.null_a`). Null B is `ccc_nulls.composite` only.
- The one CellPhoneDB-specific piece is `side_score`, about 6 lines built from `cpdb_null.feature_stats`
  and `cpdb_null.pair_scores`. Nothing in `tools/cellphonedb` is changed.
- Both nulls use the same rows, scores (checked identical on every run) and permutations.

```bash
python -c "from cellphonedb.utils import db_utils; db_utils.download_database('data/cellphonedb_v5.0.0', 'v5.0.0')"
DB=data/cellphonedb_v5.0.0/cellphonedb.zip
python -m benchmarks.null_b_cellphonedb --db $DB --part calibration --tissue mammary   # and --tissue liver
python -m benchmarks.null_b_cellphonedb --db $DB --part sensitivity --tissue mammary
python -m benchmarks.null_b_cellphonedb --db $DB --part tnbc
```

## Calibration (Tabula Sapiens, all cells): known-null data

Receptor-side (or ligand-side) genes are moved to random cells by one permutation. Every pair whose
receptor (ligand) partner is entirely shuffled is null for "s signals to t". There are 20 replicates
per scenario. FDR is the mean of V / R over replicates.

| FDR at α = 0.05 (calls per replicate) | Mammary, receptor | Mammary, ligand | Liver, receptor | Liver, ligand |
|---|---|---|---|---|
| Null A | 1.00 (6,509) | 0.97 (6,597) | 0.99 (18,433) | 0.96 (16,687) |
| **Null B** | **0 (0)** | **0 (0)** | **0 (0)** | **0 (0)** |

- **Type-I error on known-null rows (replicate 0):**

  | | p ≤ 0.01 | p ≤ 0.05 | p ≤ 0.10 |
  |---|---|---|---|
  | Null A | 36–38% | 42–45% | 46–49% |
  | Null B | 0.5–0.9% | 2.6–3.6% | 6.0–7.3% |

  Null B's p-values are valid and slightly conservative.
- **α = 0.01 and 0.10, and Benjamini–Yekutieli:** same pattern.
- **Fully shuffled labels (10 replicates per tissue):** 0 calls from both nulls (P(any call) = 0, CI 0–0.31). This is a sanity check only; it cannot tell the nulls apart.
- **B = 500, 1,000, 5,000** (Mammary, receptor, replicate 0): Null A makes 6,538, 6,520 and 6,586 calls, all false; Null B makes 0 each time.
- **Caveat:** in CellPhoneDB almost every receptor gene is only ever a receptor (544 of 549 receptor partners in Mammary). Shuffling one side therefore makes nearly every row a known null, an easier test than CellChat's (Null B FDR 0–3.1% there).

## Real labels and TNBC (held out, run once)

| | Mammary | Liver | TNBC |
|---|---|---|---|
| Null A calls | 6,209 | 15,768 | 14,780 |
| Null B calls | 2,902 (47%) | 7,170 (45%) | 9,518 (64%) |
| Null B calls outside Null A's | 0 | 0 | 0 |

- **Calls only Null A makes:**
  - The side that fails Null B splits evenly between ligand and receptor (TNBC: 2,590 vs 2,672).
  - They are spread over all ligand-expression deciles, so this is not an abundance effect, unlike CellChat (69% in the top decile).
- **Open Problems cytokine truth (TNBC, 1,000 bootstrap resamples):**
  - AUPRC: Null B 0.242, Null A 0.278 (difference −0.036, 95% CI −0.078 to 0.003), CellPhoneDB mean ungated 0.229, random 0.214.
  - Odds ratio: Null B 2.20, Null A 4.29.
  - No gain.
- **Fake gene pairs (TNBC):** 910 stratum-matched non-CellPhoneDB gene pairs. Null B cuts fake calls from 6,033 to 1,349 (−78%); fakes go from 19.4% to 7.6% of calls.
- **Patient-half reproducibility:** Jaccard 0.34 for Null B vs 0.41 for Null A. Null B is less reproducible, as on CellChat.
- **Not usable here:** `null_b|testable` (BH over every row that some labelling scores) makes 0 calls at α = 0.05. With B = 1,000 the p-value floor (1/1,001) is too coarse for a family of about 360,000 rows. Its AUPRC of 0.602 is an artefact of having no calls.

## Against the shared checklist

| Points | Status |
|---|---|
| A1, U: hypothesis, statistical claim first | met |
| B2–B5: one side permuted, sizes and expression fixed, no strata | met |
| C6, D7, D8: Null A keeps counts; same candidates, score and permutations | met (checked on every run) |
| E9: all tool adapters | partly: CellChat and CellPhoneDB on one shared Null B; NATMI and Connectome not yet |
| E10, F11: common schema, exact tool score | met (`combine`; CellPhoneDB's own means) |
| G12: number of permutations | met (500, 1,000, 5,000) |
| H13, I14: fully permuted data, empirical FDR with CIs | met |
| J15: q-value procedure | met (BH over the tested rows; BY as a check) |
| K, L: primary endpoint, A-but-not-B analysis | met |
| S: leakage | partly: the expression gate uses the observed labels, but calibration shows no inflation |
| T: sensitivity | partly: B and α done; δ and HDMT not run for CellPhoneDB |
| Level 1: calibration | met on both tissues |
| Level 2: biological check | partly: cytokine truth only, no gain |
| M, N, Level 3: NicheNet / Cauchy | shared pipeline |

## A note for the CellPhoneDB adapter

`combine(check_scores=True)` allows 5.01e-4 between the recomputed score and CellPhoneDB's 3-decimal
mean. CellPhoneDB averages in float32, so values on a rounding boundary differ by up to 5.09e-4 here.
1–3 rows exceeded the tolerance in 29 of the 103 runs, and the check stops the run. The script runs the
same check with a float32-aware tolerance, 5e-4 + 2e-5 × |score| (every run passes).
