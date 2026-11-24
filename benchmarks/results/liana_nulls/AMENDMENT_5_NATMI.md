# NATMI through benchmark 1 (AMENDMENT_5), and the additions run with it

Plan: `benchmarks/liana_nulls/AMENDMENT_5.md`. Everything marked **addition** was not pre-registered: it was decided
after AMENDMENT_5's results were known and is reported beside them, never in place of them. Result files:
`amendment5_*.json` (pre-registered), `addendum_natmi_b5000_*.json`, `addendum_nichenet.json`,
`truth_tnbc_wu2021_xenium.csv` and `rerun/` (additions). Runs: lab server, 13 November 2026.

## How this run differs from the team's machines

- **TNBC and Allen brain** were rebuilt from the figshare files Open Problems' v1 loaders download (TNBC md5
  `48365dc8…`, as on figshare), filtered and normalised as Open Problems does (all-zero genes and cells dropped,
  ln(1 + CP10k)). Shapes match the published files (42,512 × 28,200, 29 types; 14,249 × 34,617, 23 types); the file
  md5s do not (rebuilt, not copied). Every reference file (CellChatDB, CellPhoneDB v5.0.0, CytoSig, PROGENy, HCOP)
  matches `aux_files.md5`; MERFISH matches AMENDMENT_1's SHA-256.
- **Hypothalamus** rebuilt from GEO GSE113576 + MERFISH by `hypothalamus.build()`: 26,604 cells in 8 classes
  (the team's file: 29,917). Label transfer (HVG, PCA, k-means) differs between library versions.
- **Native NATMI precision fix** (`tools/natmi/adapter.py`): `em.csv` was written with 9 digits for float32 data;
  NATMI parses it as float64, so scores differed by up to 5e-9 (relative), above the plan's 1e-9. Written with 17
  digits, the largest difference is 1e-12.
- **CellChat and CellPhoneDB** row tables (not in git) were rerun here with `rerun_tools.py` (column 1 not rerun;
  column 2 equalled it on every row in the team's T0). Calls reproduce the team's exactly on TNBC (CellChat
  13,024 / 12,972 / 8,837; CellPhoneDB 14,881 / 14,765 / 9,528).
- Permutations run on 8 processes (`CCC_NULLS_WORKERS`); p-values are identical to one process (checked on
  every row of Mammary).

## T0 fidelity and calls (B = 1,000)

| Data set | Cells | Types | Native = Python (rows; max rel. diff) | Null A (BH) | Null B (BH) | B2 | LOO-t |
|---|---|---|---|---|---|---|---|
| TNBC | 42,512 | 29 | 891,945 rows; 1.0e-12 | 116,784 | 36,486 | 73 | 31,862 |
| Allen brain | 14,249 | 23 | 533,189; 9.9e-13 | 59,331 | 21,777 | 1 | 13,905 |
| TS Mammary | 30,934 | 16 | 282,934; 1.0e-12 | 39,459 | **0** | 4 | 12,262 |
| TS Liver | 22,214 | 24 | 714,335; 1.0e-12 | 87,805 | **0** | 5 | 22,589 |
| Hypothalamus | 26,604 | 8 | 77,424; 9.6e-13 | 11,755 | 3,945 | 0 | 3,112 |

Fidelity passed on all five. LIANA+'s own NATMI agrees with ours (Spearman 0.99999999999 on 46,716–120,160 shared
rows; mouse data sets not comparable, gene names differ). Null B's zero on Mammary and Liver is the p-value floor:
NATMI averages over all cells, so its BH family is huge (555,000–1,281,000 rows) and 1/1001 is too coarse.
**Addition, B = 5,000:** Null B makes 10,369 (Mammary) and 24,929 (Liver) calls; TNBC 37,657, Allen 22,897,
hypothalamus 4,175.

## T1 calibration (known nulls; FDR at BH 0.05, mean of 20 replicates; shuffled labels: replicates with any call)

| Null | Mammary rec. / lig. | Liver rec. / lig. | Labels shuffled |
|---|---|---|---|
| Null A | 0.84 / 0.90 | 0.83 / 0.91 | 0 / 10 |
| Null B | 0 / 0 (**no calls**) | 0 / 0 (**no calls**) | 0 / 10 |
| B2 | 0 / 0 (0–3 calls) | 0 / 0 | 0 / 10 |
| LOO-t | 0.88 / 0.91 | 0.90 / 0.92 | **10 / 10** |

**Addition, B = 5,000:** the same to two decimals; Null B still makes no call. With one side's genes shuffled almost
no row stays specific on both sides, so Null B's zero is consistent with validity but cannot measure its FDR.

## T2 power (planted pairs on shuffled Mammary; 21.4 of 30 planted pairs in NATMI's database)

| Rule | f = 1.5 | f = 2 | f = 4 |
|---|---|---|---|
| Null A, BH | 0 calls | 0 | 0 |
| Null B, BH | 0 calls | 0 | 0 |
| Null A, p < 0.05 | power 1.00, FDR 0.997 (13,344 calls) | 1.00 / 0.997 | 1.00 / 0.997 |
| **Null B, p < 0.05** | **1.00, FDR 0.937 (627 calls, 99 one-sided)** | 1.00 / 0.937 | 1.00 / 0.936 |
| B2, BH | 0.34, FDR 0 | 0.58 / 0 | 0.74 / 0 |
| LOO-t, BH | 1.00, FDR 0.985 | 1.00 / 0.986 | 1.00 / 0.987 |

**Addition, B = 5,000:** Null A passes BH (power 1.00, FDR 0.975–0.982); Null B still makes no call (about 22 true
rows cannot fill BH's floor over ~260,000 rows).

## T3 biology (B = 1,000): Null A's calls kept vs dropped, odds ratio [95% bootstrap CI]

| | TNBC | Allen brain | Hypothalamus |
|---|---|---|---|
| Spatial: Null B | **0.68 [0.48, 0.89] ↓** | 0.93 [0.76, 1.11] | 1.06 [0.76, 1.42] |
| Spatial: LOO-t | 1.05 | 1.68 [0.99, 2.76] | 0.86 |
| CytoSig D1 / D2 receptor (Null B) | **1.60 ↑ / 2.46 ↑** | 0.94 / **1.77 ↑** | 0.99 / **2.64 ↑** |
| PROGENy D2 ligand (Null B) | **0.64 ↓** | 0.64 | **0.22 ↓** |

Null A's own calls sit at chance on every spatial truth (0.86–1.00× the colocalised share of all rows).
**Addition, B = 5,000 (Mammary, Liver):** CytoSig D1 1.37 / **1.51 ↑**, PROGENy D1 **1.62 ↑ / 1.79 ↑**, D2 receptor
**↑** in all four.

## Verdict under the plan's rule

- **Null B:** criterion 1 holds only vacuously (no call in any calibration or power setting); criterion 2 fails
  (TNBC spatial CI below 1); criterion 3 holds. **Valid (vacuously), not biologically distinguishable.** Extended
  (Amendment 1): downstream **against** (PROGENy D2 ligand below 1), although the pre-registered D2 receptor
  prediction holds on three data sets.
- **Null A:** FDR 0.83–0.91, **rejected**, as for CellChat and CellPhoneDB.
- **LOO-t:** FDR 0.88–0.92 and calls on every shuffled-label replicate, **rejected**.
- **B2:** FDR 0 everywhere, power 0.34–0.74 under BH; too few real-data calls for biology.

## Addition: a single-cell TNBC-type spatial truth (Xenium)

10x Xenium human breast cancer (Janesick et al. 2023, GEO GSE243275; ER+/PR−/HER2+, **not triple negative**, another
patient), 2 sections, 286,532 cells. TNBC's 29 labels and Xenium's 20 mapped to 10 shared classes; the
hypothalamus neighbourhood method (10-NN, 200 permutations, z ≥ 1.645): 36 of 90 class pairs colocalised. It barely
agrees with the Visium truth (2 shared colocalised label pairs of 162). Bootstrap over the 90 class pairs:

| Tool | Null A vs chance | Null B | B2 | LOO-t | Filter (step 3 below) |
|---|---|---|---|---|---|
| NATMI (B = 5,000) | 1.11× | **1.36 [1.19, 1.57] ↑** | 2.06 [0.50, 6.33] | 1.15 [0.95, 1.35] | **1.89 [1.34, 2.74] ↑** |
| CellChat | 1.00× | **1.54 [1.09, 2.15] ↑** | 1.53 [0.06, 5.65] | 1.34 [0.81, 2.24] | 1.56 [0.71, 3.09] |
| CellPhoneDB | 1.04× | **1.33 [1.04, 1.70] ↑** | 5.23 [0.19, 29.4] | – | 1.32 [0.71, 2.56] |

## Addition: NicheNet (v2, Zenodo 7074291) for NATMI, CellChat, CellPhoneDB

**Step 1, NicheNet as a third downstream truth** (Null B's D1; footprint = a ligand's top 250 targets):

| Tool | TNBC | Allen brain | Mammary | Liver | Hypothalamus |
|---|---|---|---|---|---|
| NATMI (B = 5,000) | 1.07 | **0.69 ↓** | **1.23 ↑** | **1.61 ↑** | 0.84 |
| CellChat | **2.08 ↑** | **1.68 ↑** | 0.91 | 0.96 | 0.36 (wide) |
| CellPhoneDB | **1.29 ↑** | **1.36 ↑** | 0.93 | 1.19 | 1.61 (wide) |

**Step 2, footprint p** (rank of the ligand among ~1,225 NicheNet ligands in the receiver, after z-scoring each ligand
across receivers): on shuffled labels 4–9% of ligands are significant in more receivers than Binomial(K, 0.05)'s 99%
bound allows (expected ≤ 1%): some ligands are favoured whatever the labels. On real calls, 5–10% have p < 0.05.

**Step 3, filter = Null B's calls whose ligand has a footprint (p < 0.05) in the receiver**, spatial truths only:

| Truth | NATMI: Null B → filter | CellChat | CellPhoneDB |
|---|---|---|---|
| TNBC Visium | 0.68 → **0.47 ↓** | 1.05 → 0.73 | 0.80 → 0.33 |
| TNBC Xenium (class-pair CI) | 1.36 → **1.89 ↑** | 1.54 → 1.56 | 1.33 → 1.32 |
| Allen brain | 0.92 → 0.85 | 1.16 → 0.70 | 1.10 → 1.29 |
| Hypothalamus MERFISH | 1.05 → 1.48 | **1.95 ↑** → 0.69 | 0.60 → 1.35 |

The filter keeps 5–10% of Null B's calls. The intersection–union form, BH on max(p_B, p_footprint), makes **no**
call anywhere (the footprint p's floor is 1/n_ligands, reached by one ligand per receiver). The Cauchy combination
(the brief's suggestion) keeps 75–95% of Null B's calls and tracks Null B everywhere: as an OR it barely filters.

## Reading

- NATMI joins CellChat and CellPhoneDB: Null A is far from calibrated; Null B is valid but, on NATMI's very large
  family, needs more than 5,000 permutations, or a smaller family, for BH to call sparse signals.
- On **single-cell** spatial truths Null B keeps the colocalised pairs (Xenium: all three tools; MERFISH: CellChat);
  on the **spot-based** Visium truths no null beats chance. This is a hypothesis from a handful of tests, and the
  Xenium tumour is not triple negative.
- The NicheNet footprint is weak evidence here: as a truth it leans towards Null B (6 ↑, 1 ↓ of 15), as a filter it
  is not supported.
