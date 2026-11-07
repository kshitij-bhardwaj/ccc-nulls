# Method

## The question a call answers

A ligand–receptor call "cell type $s$ signals to cell type $t$ through $L$–$R$" makes two claims: the ligand
side is specific to the sender $s$, **and** the receptor side is specific to the receiver $t$.

CellPhoneDB and CellChat, and LIANA+'s consensus that uses them, test a different null (**Null A**): they permute
the cell-type labels of every cell at once and recompute the score $S(L_s, R_t)$,

$$p_A = \frac{\#\{b : S(L_s^{\pi_b}, R_t^{\pi_b}) \ge S(L_s, R_t)\}}{B}.$$

Null A is rejected as soon as **either** side is specific. A ligand made only by $s$, paired with a receptor made
by every cell type, is called for every receiver $t$: a one-sided interaction.

## Null B: the composite null

Null B tests the union null of the two-part claim,

$$H_0 = H_0^L \cup H_0^R,$$

where $H_0^L$ says the sender side is unrelated to the labels (receiver side as observed) and $H_0^R$ the other
way round. Each side gets its own permutation p-value, recomputing the **tool's own score** with only that side
permuted:

$$p_L = \frac{1 + \#\{b : S(L_s^{\pi_b}, R_t) \ge S(L_s, R_t)\}}{B + 1}, \qquad
p_R = \frac{1 + \#\{b : S(L_s, R_t^{\pi_b}) \ge S(L_s, R_t)\}}{B + 1}, \qquad p_B = \max(p_L, p_R).$$

This is the intersection–union test (Berger 1982; Berger & Hsu 1996): $p_B$ is a valid level-$\alpha$ test of
$H_0$ with no multiplicity adjustment. It is the partial-conjunction p-value for "both components non-null"
(Benjamini & Heller 2008), and the joint-significance test of mediation (MacKinnon et al. 2002). A score of zero
gets $p = 1$.

## Inside LIANA+

LIANA+ scores a permutation method with `fun(x, perm_stats)`, where `perm_stats` holds the ligand and receptor
statistics under each of its label permutations. `ccc_nulls.liana_consensus` wraps that function: LIANA+'s output
is returned untouched, and $p_L$, $p_R$ are computed from the **same** permutations, so Null B costs no extra
permutations. The magnitude (`lr_means`, `lr_probs`) is always LIANA+'s.

| Tool | Score $S$ | Side statistic |
|---|---|---|
| CellPhoneDB | $(L + R)/2$, zero if either is zero | mean |
| CellChat | $LR / (K_h + LR)$ | trimean, scaled by the matrix maximum |

The consensus is then rebuilt by LIANA+'s own `_aggregate` (robust rank aggregation by default):

| Variant | What changes |
|---|---|
| `liana6` | LIANA+'s five methods plus CellChat: the baseline |
| `swap` (recommended) | CellPhoneDB's and CellChat's p replaced by $p_B$ |
| `extra` | $p_B$ added as two further methods (7th and 8th), next to LIANA+'s own p |

`swap` is recommended: in the benchmark it beats `extra` on every ranking measure, because in `extra` the tools'
own Null A p-values stay in the consensus.

## Donor structure

Permuting labels across all cells treats cells of one donor as independent. When cell types follow donors,
donor differences then pass as cell-type specificity. On labels shuffled within donor (a known null), Null B with
global permutations gave p ≤ 0.05 for 11–18% of rows on TNBC and Tabula Sapiens Liver; with permutations within
donor, 0.3% (AMENDMENT_2, A2.3).

So `rank_aggregate` permutes **within donor by default** (`strata="auto"`), for Null A and Null B alike
(`ccc_nulls.liana_consensus.stratified_permutations`):
- the donor column is the first of `donor_id`, `donor`, `patient_id`, `patient`, `individual`, `subject`,
  `orig.ident`, `sample_id`, `sample` in `adata.obs` with two or more values;
- if fewer than two donors hold two or more cell types, no label could move: global permutations, with a warning;
- without a donor column: global permutations, as LIANA+;
- `strata="<column>"` names the column; `strata=None` asks for global permutations.

The column used is in `res.attrs["strata"]`. The pre-registered benchmark runs used global permutations, as
LIANA+ does, apart from A2.3.

## What is inherited from LIANA+

For a complex, `perm_stats` follows the one subunit that is lowest under the observed labels, as in LIANA+'s own
permutation test. The fidelity tests check that the five-method consensus reproduces `li.mt.rank_aggregate`
exactly, and that the captured tools reproduce `li.mt.cellphonedb` and `li.mt.cellchat` exactly.

## Other nulls in the package

The package also contains a gene-stratum decoy null and **B2** (decoy genes with a CAMERA-style correction for
inter-gene correlation; `ccc_nulls.competitive`), which the benchmark reports for comparison, and the original
significance layer (`ccc_nulls.significance_layer`) used in the project's earlier stages.

## References

- Berger RL. Multiparameter hypothesis testing and acceptance sampling. *Technometrics* 24, 295–300 (1982).
- Berger RL, Hsu JC. Bioequivalence trials, intersection–union tests and equivalence confidence sets. *Statist. Sci.* 11, 283–319 (1996).
- Benjamini Y, Heller R. Screening for partial conjunction hypotheses. *Biometrics* 64, 1215–1222 (2008).
- MacKinnon DP et al. A comparison of methods to test mediation and other intervening variable effects. *Psychol. Methods* 7, 83–104 (2002).
- Dimitrov D et al. LIANA+ provides an all-in-one framework for cell–cell communication inference. *Nat. Cell Biol.* 26, 1613–1622 (2024).
- Dimitrov D et al. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. *Nat. Commun.* 13, 3224 (2022).
- Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Res.* 40, e133 (2012).
