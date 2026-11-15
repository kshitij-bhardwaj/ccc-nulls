# Final pipeline: our nulls inside LIANA+'s consensus — plan

Written before any evaluation of this pipeline. One check ran before it: the engine
(`ccc_nulls.liana_consensus`) on 3,000 hypothalamus cells with 100 permutations, which compared LIANA+'s
own outputs with ours (all identical). It also printed call counts at specificity_rank ≤ 0.05:
liana5 2,127; liana6 1,472; swap_null_b 1,592; extra_null_b 1,597; swap_b2 2,075; extra_b2 2,084;
extra_null_a 1,826 (10,156 rows). Those counts are disclosed here. The decision rule below does not use call counts.

## Design (decided with the project lead, 2026-11-03)

- LIANA+ 1.8.0's consensus (`rank_aggregate`, RRA) over its five methods, plus **CellChat as the sixth**.
- For the two methods with a permutation null, CellPhoneDB and CellChat, **only the p-value is replaced**.
  The magnitude (`lr_means`, `lr_probs`) stays LIANA+'s. Our nulls run on LIANA+'s own scores and resource:
  - `null_a`: LIANA+'s own p (all labels permuted together);
  - `null_b`: the composite null on LIANA+'s permutations, p = max(p_L, p_R);
  - `b2`: decoy genes with the CAMERA correction.
- Variants, each a full LIANA+ consensus:
  - `liana5`: LIANA+ as shipped;
  - `liana6`: + CellChat, the **baseline**;
  - `swap_<null>`: liana6 with the two tools' p replaced;
  - `extra_<null>`: liana6 plus CellPhoneDB·<null> and CellChat·<null> as the 7th and 8th methods.
- NATMI and Connectome have no permutation null in LIANA+; they enter unchanged. Our nulls for them come later.
- Settings: `n_perms` = 1,000, seed 1337, expr_prop 0.1, LIANA+'s consensus resource (mouseconsensus for
  mouse), data X = ln(1 + CP10k) as in `benchmarks/liana_nulls`.

## Calls and rankings

- A **call** is `specificity_rank ≤ 0.05`.
- The **ranking** is `specificity_rank`, ascending.

## Stages

- **F0 fidelity (gate).**
  - On every row of all five data sets, `liana5` must equal `li.mt.rank_aggregate` with the same settings,
    and the captured CellPhoneDB and CellChat outputs must equal LIANA+'s standalone methods on the hypothalamus.
  - If F0 fails, nothing else is reported.
- **F1 calibration.** Tabula Sapiens Mammary and Liver, 3 replicates each.
  - Scenarios:
    - genes used only as receptors shuffled across cells (one permutation for all of them);
    - genes used only as ligands, likewise;
    - labels shuffled.
  - FDR = calls whose shuffled side has every subunit shuffled ÷ calls. Shuffled labels: calls.
- **F2 power.** Planted edges, the `benchmarks/liana_nulls/power.py` design (Tabula Sapiens Mammary,
  labels shuffled, 30 pairs), at f = 1.5 and 2.0 with 5 replicates. Measures:
  - power on the planted rows;
  - FDP over non-ambiguous rows;
  - one-sided false calls.
- **F3 biology.** Among liana6's calls with s ≠ t, compare the calls the variant keeps with the calls it drops.
  - Spatial (TNBC, Allen brain, hypothalamus): odds ratio of co-location, 95% bootstrap CI over
    (s, t) pairs, plus pair-level AUPRC.
  - Downstream (all five data sets; CytoSig, PROGENy): D1 odds ratio of the receiver responding.

## Decision rule

Applied to `swap_null_b` and `extra_null_b` separately:

1. **Calibration:** mean FDR is below liana6's in all four tissue × side scenarios.
2. **Power:** at f = 2.0, power ≥ 0.8.
3. **Biology:** the spatial odds-ratio CI lies above 1 in at least one data set and below 1 in none.

Verdict:
- **"improves LIANA+"** if 1–3 all hold;
- **"valid, not biologically distinguishable"** if 1 and 2 hold without 3;
- **"not supported"** otherwise.

The `null_a` and `b2` variants are reported, not judged.

The replicate counts are smaller than in `benchmarks/liana_nulls` because every replicate runs all six
methods with 1,000 permutations for two of them.
