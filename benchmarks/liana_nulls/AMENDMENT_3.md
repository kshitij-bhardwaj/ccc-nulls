# Amendment 3: small effects and more permutations; two more Null B columns; an ablation

Added on 1 November 2026, after `PLAN.md`, `AMENDMENT_1.md` and `AMENDMENT_2.md`. This is an
**addition**. It changes neither the pre-registered verdict nor the extended verdict of Amendment 1.
Result files written by these parts record the SHA-256 of this file.

## What had been seen when this was written

Everything in `results/liana_nulls/README.md`: T0 to T3, Amendments 1 and 2. In particular:
- in T2, power is flat from f = 1.5 to 4;
- BH Null B makes almost no calls (CellChat about 10, CellPhoneDB 0).

Nothing below had been computed.

## A3.1 T2 at small effects and with more permutations

The T2 design of `PLAN.md` is unchanged: same 5 replicates, backgrounds, planted pairs and seeds.
Added:
- **Small effects:** f = 1.1, 1.2 and 1.3 at B = 1,000.
- **Check:** f = 1.5 at B = 1,000 is run again. The Null A and Null B columns must equal those in
  `power.json`.
- **More permutations:** B = 5,000 at f = 1.2 and 1.5.
  - Reason: with B = 1,000 the smallest p is 1/1001, so BH can reject only when about m/50 rows sit
    at that floor (Amendment 2).
  - At B = 5,000 the floor is 1/5001.

Every variant of A3.4 is scored on every run. Reported as before: power, FDR and one-sided false
calls, each a mean with a t-CI over the replicates.

## A3.2 Column "Null B (#9)": CellPhoneDB's Null B as pushed in pull request #9

`benchmarks/null_b_cellphonedb.run` is run unchanged:
- its side score: float64 means from `cpdb_null.feature_stats`;
- its permutations: `default_rng(0)`, B = 1,000, ties counted;
- its three procedures: BH over CellPhoneDB's gate (primary), BH over `composite.testable`, and
  Benjamini–Yekutieli (BY).

CellPhoneDB only, since #9 has no CellChat path. For mouse data the genes are first renamed to
human orthologs, as column 3 does. Reported against our column 3:
- the share of rows with the same call;
- the Jaccard index of the call sets;
- the Spearman correlation of p.

## A3.3 Column "Null B2": decoy genes with a CAMERA-style correction

**Why.**
- The October decoy null (`ccc_nulls.nulls.null_b`) ranked a candidate against abundance-matched
  non-database genes. On shuffled labels it gave p < 0.05 for 13–14% of tests.
- Its premise is that the decoys are independent draws from the candidate's null. That is false:
  genes are correlated across cells.
- This is the known flaw of competitive gene-set tests by gene sampling (Goeman & Bühlmann 2007).
- CAMERA (Wu & Smyth 2012, *Nucleic Acids Res.* 40:e133) keeps the competitive null and puts the
  inter-gene correlation into the variance. Efron & Tibshirani's restandardization (2007) is the
  same idea.

**Null hypothesis of one side.** The side's genes are no more specific to the cell type than
comparable genes that are not signalling genes.

**Definition.**
1. **Gene statistic:** the one-vs-rest Welch t of gene g in type k on ln(1 + CP10k), for all genes.
2. **Side of a row:**
   - the ligand side is the ligand's subunit genes S in the source;
   - the receptor side is the receptor's subunit genes in the target;
   - CellChat's agonists, antagonists and co-receptors are not included.
3. **Gate:** every subunit is detected in ≥ 10% of the type's cells. Otherwise the side's p is 1.
4. **Decoys:**
   - genes in no database used here (CellChatDB, CellPhoneDB, LIANA+ consensus of the species);
   - in the subunit's abundance stratum (`gene_strata`: 10 mean bins × 3 detection bins over the
     loaded genes);
   - detected in ≥ 10% of the type's cells, the same gate as the candidate;
   - at most 200 per (stratum, type), drawn with seed 0.
   - The reference D of a side is the union over its subunits. For Tabula Sapiens, where only
     database genes were loaded, 5,000 non-database genes drawn with seed 0 are loaded as well.
     Library sizes, and so values, are unchanged.
5. **Correlation:** Pearson correlation of within-type residuals (a cell's value minus its type's
   mean), on 5,000 cells drawn with seed 0 (all cells if fewer). This is CAMERA's residual
   correlation.
6. **Statistic:** T = mean over S of t minus mean over D of t.
7. **Variance under the competitive null:**
   - Var(T) = τ² (VIF_S/m + 1/n) + VIF_S/m + VIF_D/n − 2 r̄_SD, where m = |S| and n = |D|;
   - VIF_S = 1 + (m − 1) ρ̄_S, CAMERA's variance inflation factor (ρ̄_S is the mean correlation
     among S);
   - VIF_D = 1 + (n − 1) ρ̄_D (ρ̄_D is the mean correlation among the decoys);
   - r̄_SD is the mean correlation between S and D;
   - τ² = max(0, s_D² − (1 − ρ̄_D)), where s_D² is the variance of the decoys' t across genes. It
     is the spread of real specificity among comparable genes, after the part that correlation
     explains is removed.
8. **p-value:**
   - side p is the upper tail of t with n − 1 degrees of freedom at T / √Var(T);
   - row p = max(p_ligand, p_receptor), the intersection–union rule of Null B;
   - BH over the tool's family, the same rows as Null A and Null B.

**Validity check, fixed now.** On Tabula Sapiens Mammary and Liver, 10 replicates each:
- **Labels shuffled:**
  - side-level type-I error at α = 0.05 (unique gated sides);
  - P(any BH call), with a Clopper–Pearson CI.
- **Receptor genes shuffled across cells** (genes that are receptor subunits and never ligand
  subunits, one permutation): FDR = mean V/R, where false rows have a receptor side made only of
  shuffled genes.

B2 is called **calibrated** if, in every check:
- the 95% t-CI of the side-level type-I error reaches down to 0.05 or below;
- the FDR's upper CI bound is ≤ 0.05.

## A3.4 Ablation, on all five data sets

Each variant changes one part of the full null and keeps the rest.

| Null | Part | Full | Variants |
|---|---|---|---|
| Null B | sides | both, p = max | ligand side only; receptor side only; joint shuffle (= Null A) |
| Null B | multiplicity | BH | BY; none (the tool's own rule) |
| Null B | family (CellPhoneDB) | its gate | `composite.testable` |
| Null B | implementation (CellPhoneDB) | ours (its permutations, float32) | #9 (A3.2) |
| Null B | permutations | B = 1,000 | B = 5,000 (T2 only) |
| B2 | correction | CAMERA | none: Var(T) = s_D² (1/m + 1/n), the old premise |
| B2 | decoys | abundance-matched | any non-database gene passing the gate (cap 200 per type) |
| B2 | combination | both sides, p = max | joint: z = (z_L + z_R)/√2, one-sided |

Endpoints, per variant:
- **T2 (A3.1 runs):** power, FDR, one-sided false calls.
- **Shuffled labels and receptor-shuffled data (B2 variants):** as in the A3.3 validity check.
- **Real data, each data set:**
  - calls;
  - the Jaccard index with the full Null B;
  - the calls the variant makes that Null A does not.
- **Biology:**
  - the kept-vs-dropped odds ratio, as in `PLAN.md`, where "kept" means a Null A call that the
    variant also calls;
  - spatial truth on TNBC, Allen brain and the hypothalamus;
  - downstream D1 (CytoSig, PROGENy) on all five data sets;
  - the same bootstraps as before.

The ablation has no decision rule. Two things are stated in advance:
- if the ligand-side-only variant's T2 FDR is close to Null A's, the receptor side is what controls
  the false calls (and the reverse);
- if B2 without the correction fails its validity check and B2 with it passes, the correction is
  what repairs the decoy null.

## Pipeline

- B2 lives in the shared package: `ccc_nulls.competitive`.
- `core.run_dataset` writes `p/q/call` columns for B2 (and, for CellPhoneDB, for #9) into
  `liana_res`.
- The part `extend` adds them to cached core results without recomputing columns 1–3.

## Status

The report's status table gains one row per A3 component.
