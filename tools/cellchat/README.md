# CellChat

One owner, assigned within the team. Only the owner edits this folder; anything shared goes in
`src/ccc_nulls/` through a reviewed pull request.

| | |
|---|---|
| Tool | CellChat **2.1.2**, the latest release (`jinworks/CellChat`, tag `v2.1.2`, April 2025), **ported to Python**: no R needed |
| Also available | `2.2.0.9001` (GitHub main, `75253cd`) and `1.0.0` (the 2021 paper, `sqjin/CellChat`, tag `v1.0.0`) |
| Statistic tested by the nulls | `prob`, CellChat's communication probability (LIANA's re-implementation calls it `lr_probs`) |
| Native null (Null A) | Yes: `computeCommunProb` permutes the cell labels, `nboot = 100`, `seed.use = 1` |
| Extra Python dependency | `rdata` (reads CellChatDB from the package's `.rda` once): `pip install rdata` |

## Files

| File | What it is |
|---|---|
| `cellchat_py.py` | CellChat ported line by line from R: CellChatDB, `subsetDB`, `subsetData`, `identifyOverExpressedGenes` (presto's Wilcoxon test), `identifyOverExpressedInteractions`, `computeCommunProb` (2.1.2, 2.2.0.9001, 1.0.0), and R's random-number generator (`set.seed`, `sample.int`) |
| `adapter.py` | dataset → CellChat inputs (`prepare_inputs`) → CellChat run (`run_native`, the port) → common schema (`native_to_schema`) |
| `nulls.py` | `null_a`: CellChat's permutation test, vectorised; `null_b`: the composite null (max of the sender-side and receiver-side permutation p); `testable`: the multiple-testing family; `null_b_empirical`: B1, kept for the negative result |
| `match_null_a.py` | TNBC: port vs Null A for every version, independent permutations, version differences, ground-truth metrics with bootstrap CIs |
| `calibrate_null_a.py` | Tabula Sapiens (Mammary, Liver): full-permutation calibration |
| `null_b_calibration.py` | Tabula Sapiens: Null B vs Null A on known-null data (labels, receptor-side and ligand-side gene shuffles) and on the real labels |
| `null_b_tnbc.py` | TNBC, held out: Null B vs Null A against the cytokine truth, disagreements, patient-half reproducibility |
| `tail_pvalues.py`, `empirical_null.py` | B1: randomized permutation p, GPD and Pearson III tail p-values, Efron's empirical null |
| `b1_calibration.py`, `b1_tail_accuracy.py`, `b1_diagnosis.py` | B1: calibration, tail accuracy against 100,000 shuffles, why the empirical null is not identifiable |
| `results/` | the JSON outputs of the scripts |

## Versions

All three use the same CellChatDB file (identical in v2.1.2 and 75253cd) and the same filters.

| | 1.0.0 (paper) | **2.1.2 (release, default)** | 2.2.0.9001 (main) |
|---|---|---|---|
| complexes, co-receptors | per cell, then trimean per group | trimean per gene, then combined | as 2.1.2 |
| quantile in `triMean` | `stats::quantile`: `(1 - h) a + h b`, only where `b != a` | `stats::quantile` | `collapse::fquantile`: `a + h (b - a)` |
| p-value | `#(prob_perm >= prob) / nboot` | `#(prob_perm > prob) / nboot` | as 2.1.2 |
| default data | smoothed over a protein network (`projectData`); **we run it on raw data** | raw | raw |

2.1.2 and 2.2.0.9001 differ only in the last bit of some quantiles: on TNBC their probabilities
differ in 1,186 of 1,757,690 rows by at most 6 × 10⁻¹⁷, and every p-value is the same.

## What CellChat computes, and what Null A re-implements

From `computeCommunProb` (R/modeling.R, v2.1.2), on `data = data.signaling / max(data.signaling)`:

| Step | CellChat | ours |
|---|---|---|
| group average | `triMean`: `mean(stats::quantile(x, c(.25, .5, .5, .75)))` per gene and group | same order statistics and interpolation, R's two-pass `mean` |
| ligand / receptor | single gene: its average; complex: geometric mean of the subunits' averages | same |
| co-receptors | `R * prod(1 + avg[co-A]) / prod(1 + avg[co-I])` | same, sequential products |
| law of mass action + Hill | `P1 = LR / (Kh + LR)`, `LR[s, t] = L[s] R[t]`, `Kh = 0.5`, `n = 1` | same |
| agonist / antagonist | `P2 = ag agᵀ`, `ag = prod(1 + avg / (Kh + avg))`; `P3 = an anᵀ`, `an = prod(Kh / (Kh + avg))` | same |
| population size | `P4`, off by default (`population.size = FALSE`) | `population_size=False` |
| **Null A** | `set.seed(seed.use)`, `nboot` permutations of all cells' labels (`sample.int`), recompute `prob`; `p = #(prob_perm > prob) / nboot`; `p = 1` where `prob = 0` | same; `rng="cellchat"` draws CellChat's own permutations with R's generator |
| significance | `subsetCommunication(thresh = 0.05)`: `p < 0.05`. **No q-values, no FDR control** | BH q-values added, for the comparison |

- CellChat's p has no `+1`, so it can be exactly 0, and with `nboot = 100` its resolution is
  0.01. `null_a` returns `p_null_a_cellchat` (CellChat's formula) and `p_null_a`
  (`(1 + #) / (B + 1)`, the repository's convention).
- Before Null A, CellChat keeps genes over-expressed in at least one group (Wilcoxon p < 0.05,
  logFC > 0) and the L-R pairs whose ligand and receptor both pass.
- LIANA's `cellchat` is a simplified re-implementation (no co-receptors, agonists or
  antagonists; LIANA's resource), so it is a reference, not the original.

## How the match is checked

```bash
python -m tools.cellchat.match_null_a --quick     # 8,000 TNBC cells
python -m tools.cellchat.match_null_a             # all 42,512 TNBC cells (~15 min)
pytest tools/cellchat                             # includes exact-match tests on 700 PBMCs
```

1. **R's generator.** `RRandom` reproduces known R outputs (`set.seed(1); sample(10)` is
   `9 4 7 1 2 5 3 10 6 8`, and four more), so the permutations are the ones CellChat draws.
2. **Reference.** `cellchat_py.compute_commun_prob` follows the R loop literally (dense sort per
   group, one L-R pair and one permutation at a time).
3. **Our Null A.** `nulls.null_a(..., B=100, seed=1)` is a separate, vectorised implementation
   (sparse order statistics, all pairs at once). It must match the reference exactly.
4. **Same null distribution.** `null_a(..., rng="numpy")` with independent permutations must
   agree with CellChat up to Monte-Carlo error.

### How far the match goes

| Claim | Established? | Evidence |
|---|---|---|
| our fast Null A = our line-by-line port | **yes, exactly** | prob, p, BH q and calls identical on every TNBC row, for all three versions |
| our permutations = the ones R draws for `set.seed(1)` | **yes** | `RRandom` reproduces 5 known R outputs |
| the port = CellChat's R source | **by reading the code** | every formula and operation order follows R/modeling.R of the version named |
| our numbers = what CellChat prints when run in R | **not checked** | CellChat in R was never run (the project is Python only) |

What could still differ from running CellChat in R, and why it should not:
- `log` / `exp`: numpy matches the C library bit for bit (checked on 10⁶ values each).
- R's `mean` uses long double; on Apple Silicon long double is double, which the port assumes.
  On x86 Linux R's `mean` may differ in the last bit, which can only flip an exact tie between a
  permuted and an observed probability.
- presto's p-values at exactly the 0.05 boundary can only change which L-R pairs enter, not Null A.

## Results: Null A on TNBC (Open Problems `tnbc_wu2021`, 42,512 cells, 29 groups)

`results/null_a_match_tnbc_wu2021_full.json` (`_quick` for 8,000 cells). CellChat's pipeline:
1,147 CellChatDB genes in the data, 1,143 over-expressed, 2,090 L-R pairs, 29 × 29 × 2,090 =
1,757,690 (pair, source, target) rows. The trimean is 0 unless more than ~25% of a group
expresses the gene, so most rows are 0 and get p = 1.

| Check | 2.1.2 (release) | 2.2.0.9001 | 1.0.0 |
|---|---|---|---|
| rows with prob > 0 | 14,872 | 14,872 | 14,031 |
| prob identical, port vs Null A (all rows) | **1,757,690 / 1,757,690** | 1,757,690 / 1,757,690 | 1,757,690 / 1,757,690 |
| p identical (rows with prob > 0) | **14,872 / 14,872** | 14,872 / 14,872 | 14,031 / 14,031 |
| BH q identical | **14,872 / 14,872** | 14,872 / 14,872 | 14,031 / 14,031 |
| calls at p < 0.05, ours / CellChat's | 13,024 / 13,024, 0 differ | 13,024, 0 differ | 12,470, 0 differ |
| p = 0 (no permutation reached the observed prob) | 12,862 | 12,862 | 12,300 |

Release 2.1.2 with independent permutations (same null, different draws):

| | within 2 Monte-Carlo sd | within 3 sd | calls at p < 0.05 (both / only ours / only CellChat) |
|---|---|---|---|
| B = 100 | 99.8% | 99.99% | 13,002 / 17 / 22 |
| B = 1,000 | 99.6% | 100% | 13,007 / 18 / 17 |

1.0.0 vs 2.1.2: 3,260 of the 1,757,690 probabilities differ, all involving complexes or
co-receptors; 13,765 of the 14,031 shared p-values are identical and 13,830 are ≥ 2.1.2's (ties
count against in 1.0.0); calls at p < 0.05: 12,371 in both, 653 only in 2.1.2, 99 only in 1.0.0.

CellChat's own significance on TNBC: 13,024 of 14,872 scored interactions (88%) have p < 0.05,
12,862 have p = 0, so their BH q-values are 0 as well and BH removes nothing. CellChat reports
no FDR, and its p-value cannot be read as one.

## Evaluation against ground truth (plan step 1, TNBC cytokine truth)

Field `openproblems` of the same file. Open Problems ligand-target protocol (812 truth rows, 166
positive); point estimates equal `ccc_nulls.evaluation.openproblems`; 95% CIs from 1,000
bootstrap resamples of the truth rows, paired across methods.

| | truth rows scored > 0 | AUPRC [95% CI] | odds ratio, top 5% [95% CI] | AUPRC, 2.1.2 minus this [95% CI] |
|---|---|---|---|---|
| **CellChat 2.1.2** | 41 | **0.350** [0.258, 0.447] | 3.45 [1.73, 6.43] | |
| CellChat 2.1.2, score kept where Null A p ≤ 0.05 | 41 | 0.350 [0.258, 0.447] | 3.45 | 0 |
| CellChat 2.2.0.9001 | 41 | 0.350 | 3.45 | 0 |
| CellChat 1.0.0 (paper) | 25 | 0.291 [0.193, 0.404] | 1.95 [0.77, 3.31] | 0.059 [0.001, 0.125] |
| LIANA's CellChat re-implementation | 158 | 0.235 [0.189, 0.293] | 3.09 [1.32, 5.54] | 0.115 [0.039, 0.191] |
| random | 812 | 0.214 [0.175, 0.264] | 1.51 [0.57, 3.08] | 0.136 [0.040, 0.231] |

- CellChat 2.1.2 ranks the cytokine truth better than LIANA's re-implementation and than random
  (the paired CIs exclude 0), and slightly better than 1.0.0 (CI just above 0).
- **Caveat:** CellChat gives a score > 0 to only 41 of the 812 truth rows (1.0.0: 25); the other
  ~770 are tied at the bottom. Its AUPRC rests on those few rows, hence the wide CI, and the top
  5% (40 rows) is essentially its whole scored set.
- Gating by Null A changes nothing: every top-scoring interaction already has p < 0.05.

## Calibration (plan step 2, Tabula Sapiens, full permutation)

`python -m tools.cellchat.calibrate_null_a --tissue all [--version 1.0.0]` →
`results/calibration_<tissue>_<version>.json`. Tabula Sapiens (the root README's calibration
collection), raw counts → ln(1 + CP10k), cell types with ≥ 5 cells:
**Mammary** 30,934 cells, 16 types, 2 donors; **Liver** 22,214 cells, 24 types, 4 donors.
Every cell's label is shuffled once (20 independent shuffles per tissue), so no interaction is
real, and CellChat runs end to end with its defaults (nboot = 100, seed.use = 1). Rates are
pooled over shuffles; 95% CIs from a bootstrap over shuffles.

P(p < 0.05) on the rows CellChat reports (prob > 0); a calibrated test gives 5/101 = **4.95%**:

| Fully shuffled labels | Mammary, 2.1.2 | Liver, 2.1.2 | Mammary, 1.0.0 | Liver, 1.0.0 |
|---|---|---|---|---|
| CellChat pipeline | **15.1%** [13.7, 16.4] | **12.6%** [10.9, 14.3] | **15.7%** [14.3, 17.2] | **12.9%** [11.0, 14.9] |
| without the over-expression filter | 9.4% [8.4, 10.3] | 9.4% [8.6, 10.3] | 9.8% [8.8, 10.9] | 9.7% [8.9, 10.6] |
| rows where no permutation reached prob = 0 (no filter, 3 shuffles) | **4.92%** | **5.10%** | | |
| P(p = 0), pipeline (expected 0.99%) | 4.3% [3.7, 4.8] | 3.2% [2.6, 3.8] | 4.5% [3.8, 5.2] | 3.2% [2.6, 4.0] |
| P(p < 0.05) over all rows (prob = 0 rows have p = 1) | 0.47% | 0.16% | 0.42% | 0.13% |
| BH false calls over 20 shuffles: CellChat's p / (1 + #)/(B + 1) | 3,667 / 0 | 4,157 / 0 | 3,283 / 0 | 3,524 / 0 |

With the real labels, 79% (Mammary) and 84% (Liver) of scored rows have p < 0.05.

What this shows:
1. **Our Null A is correct.** Where nothing is selected it is calibrated (4.92% and 5.10% vs
   4.95%), and over all rows it is conservative, as a permutation test must be.
2. **CellChat's reported p-values are anti-conservative, 2.5 to 3.2 times nominal, in both tissues
   and both versions**, because rows are selected with the same labels that are tested:
   - the over-expression filter keeps genes that look high in some group under these labels
     (15.1% → 9.4% on Mammary without it);
   - only rows with prob > 0 are reported, and a permutation that drives prob to 0 always loses
     to an observed prob > 0 (11.5% on Mammary where that happens, 4.9% where it does not).
3. **CellChat's p = 0 breaks BH.** With p = 0 allowed, BH calls 164 to 208 interactions per
   shuffled data set (mean over 20 shuffles) with no signal at all; with `(1 + #) / (B + 1)` it calls none.

Full permutation separates a valid test from an invalid one here only because CellChat's
selection is label-dependent; for abundance confounding, the fake-L-R-resource test is still
needed (docs/PLAN.md).

## B1, Efron's empirical null: a negative result

B1 was the first Null B: keep Null A's statistic, but estimate the null distribution from the bulk of
the interactions' z-values (Efron 2004; 2007; 2010) instead of taking N(0, 1). It cannot be used with
CellChat. The evidence, all on Tabula Sapiens (TNBC was never touched):

**1. The z-value has to come from the permutation p.** The standardized score (prob − permutation
mean) / permutation sd is not N(0, 1) under the null: prob is 0 in most permutations, so a rare
non-zero value looks extreme. On fully shuffled Mammary (`b1_calibration.py`, 2 shuffles, B = 1,000,
`results/b1_calibration_mammary_smoke.json`), testing it against N(0, 1) gives 271 and 72 false
calls at q ≤ 0.05, and the empirical null on it 1,534 and 824. The probit of the permutation p,
z = Φ⁻¹(1 − p) with ties broken at random (`tail_pvalues.randomized_p`, exactly U(0, 1) under the
null), gives 0 and 0 under both nulls.

**2. Tail p-values below 1/B: two methods, GPD chosen.** `tail_pvalues.py`, known-answer tests in
`test_b1_empirical.py`, rule fixed before the results: calibration first, then tail accuracy, ties
to GPD.

| | Generalized Pareto tail (Knijnenburg et al. 2009), shape ξ ≥ 0 | Pearson type III, three permutation moments (Mielke & Berry 2007) |
|---|---|---|
| synthetic tail, gamma (Pearson III's own family), p = 10⁻⁴, 10⁻⁵ | 2–4× too large (conservative) | exact |
| synthetic tail, lognormal | within ~3× | 24× too small at 10⁻⁵ |
| shuffled Mammary, false calls (2 shuffles) | 0, 0 | 10, 0 |
| vs 100,000 permutations, real Mammary labels, p_ref in [10⁻⁴, 10⁻³) (22 rows) | median 2.3× too large; 14% more than 2× too small | median exact; **41%** more than 2× too small |
| same, p_ref in [10⁻³, 10⁻²) (45 rows) | 1.3× too large; 9% | exact; 20% |

`b1_tail_accuracy.py` → `results/b1_tail_accuracy_mammary.json` (40 random L-R pairs, B = 1,000
estimates vs B = 100,000, Knijnenburg et al.'s design; 600 further tail rows lie below 10⁻⁴, beyond
the reference). Holding ξ ≥ 0 is our one change to Knijnenburg et al.: with a few hundred exceedances
the shape MLE is biased below 0 for exponential-type tails, and ξ < 0 implies an upper endpoint
past which p = 0 (median 20× too small at p = 10⁻⁴ without it).

**3. The empirical null is not identifiable for CellChat.** `b1_diagnosis.py` →
`results/b1_diagnosis_<tissue>.json`. Efron's estimate needs the central z's to be null; he calls it
reliable when about 90% or more of the tests are null (2010, ch. 6). Here:

| | Mammary | Liver |
|---|---|---|
| testable rows with prob > 0, observed / under random labels | **14.5% / 35.0%** | **23.2% / 33.5%** |
| reported rows Null A calls (BH, 0.05) | 7,095 of 8,952 (79%) | 13,314 of 15,839 (84%) |
| reported rows whose p is below what 1,000 permutations resolve (< 30 non-zero permutations, no tail fit) | 4,790 | 8,098 |
| calls at q ≤ 0.05, empirical null on all testable rows, over 6 standard settings (MLE / central matching × central 80%, 50%, 33%) | **1,477 to 7,096** | **2,101 to 13,938** |
| same, reported rows only, each compared with its own non-zero permutations | **22 to 932**; central matching fails (log-density not concave) | **325 to 2,150** |

Cell types differ in almost every gene, so most (pair, sender, receiver) rows are non-null, mostly in
the direction of absence (fewer non-zero rows than random labels give). There is no null bulk: the
reported rows' z-values run from −3.9 to +30 with median 2.4 (Mammary). The fitted "null" then
depends on the tuning setting, not the data, and the calls vary 40-fold between settings on the same
data set. B1 was therefore not frozen and not run on TNBC.

What B1 leaves behind: the probit-p construction, the tail p-values, and the lesson that Null A's
hypothesis is false almost everywhere, so a useful Null B must change the hypothesis tested, not
re-estimate the null of the same one. That is the Null B below.

## Null B: a composite null for "s signals to t" (`nulls.null_b`)

**The hypothesis.** Null A permutes all labels at once, so it tests "neither partner depends on cell
type". Rejecting it shows that the ligand *or* the receptor varies with cell type. A user who reads a
row as "s signals to t" is claiming both: the sender side is specific to s **and** the receiver side
is specific to t. Null B tests that claim directly, as the union of two nulls:

H0 = H0_L ∪ H0_R, where H0_L is "the sender side is not specific to s" and H0_R is "the receiver side
is not specific to t".

**The test (intersection-union; Berger 1982, Berger & Hsu 1996).**
- p_L: shuffle the cell labels for the sender side only (ligand, plus the sender's agonists and
  antagonists), and keep the receiver side at its observed group averages.
- p_R: the mirror image (receptor, co-receptors, and the receiver's agonists and antagonists).
- p_B = max(p_L, p_R). This is a valid p-value for the union null with no multiplicity
  adjustment, and BH over p_B controls FDR for the composite hypotheses (Benjamini & Heller 2008).
  The same construction is the "joint significance" test of mediation analysis (MacKinnon et al.
  2002); HDMT (Dai et al. 2022) and DACT (Liu et al. 2022) are power refinements of it and were not
  needed here.
- Everything else is Null A's: CellChat 2.1.2's exact score, candidates, B = 1,000 and permutation
  p = (1 + #{prob_b > prob}) / (B + 1). Both nulls use the same permutations for a given seed, and
  `score_null_b` equals `score_null_a`. The non-permuted side keeps its observed profile, and there are
  no gene bins or decoys.
- Interval null (sensitivity): p_L and p_R at δ count permutations with prob_b > prob − δ, which tests
  "specific by at least δ" (Hodges & Lehmann 1954; TREAT, McCarthy & Smyth 2009). δ ∈ {0.01, 0.05}.

**The multiple-testing family.** CellChat reports only rows with prob > 0, a set chosen on the
observed labels. BH over that set is biased: the observed labelling is the only one guaranteed to
give prob > 0, so the selection conditions on the outcome (calibration below, row "reported rows").
Null B corrects over the **testable** rows instead (`nulls.testable`). A row is testable when some
labelling gives prob > 0: the observed one, a Null A permutation, or a permutation of either side.
Membership does not depend on which labelling is observed. Rows with prob = 0 get p = 1, as in Null A.

### Calibration on known-null data (Tabula Sapiens, `null_b_calibration.py`)

Fully shuffled labels make both partners null, so they cannot separate the two tests. New known-null
data that do separate them:
- **receptor**: every receptor-side gene (receptors, co-receptors) is shuffled across cells with a
  single cell permutation. Complexes stay intact, the ligand side keeps its cell-type signal, and every
  pair whose receiver side was shuffled is null for "s → t" in every (sender, receiver).
- **ligand**: the mirror image.
- Genes on both sides (agonists, antagonists, and genes that are a ligand in one pair and a receptor
  in another) are never shuffled, and the pairs that use them are not counted as known null.

Each scenario has 20 replicates per tissue. FDR = E[V / R], where V counts calls on known-null rows
and R counts all calls. It is estimated as the mean over replicates with a t-interval, and it is a
lower bound, since calls on other rows may be false too. BH is run at α = 0.05 over the testable rows.

| FDR at α = 0.05 (95% CI) | Mammary, receptor | Mammary, ligand | Liver, receptor | Liver, ligand |
|---|---|---|---|---|
| Null A | 0.271 (0.232–0.311) | 0.717 (0.692–0.742) | 0.643 (0.621–0.664) | 0.746 (0.720–0.772) |
| **Null B** | **0.012 (0.008–0.017)** | **0 (no false call in 20)** | **0.031 (0.019–0.042)** | **0.029 (0.013–0.045)** |
| sender side only (p_L) | 0.276 | 0.024 | 0.688 | 0.041 |
| receiver side only (p_R) | 0.018 | 0.789 | 0.039 | 0.791 |
| Null A, BH over reported rows | 0.280 | 0.732 | 0.657 | 0.759 |
| Null B, BH over reported rows | 0.043 | 0.094 | 0.090 | 0.121 |
| Null B, δ = 0.01 or 0.05 | 0 (0 calls) | 0 (0 calls) | 0 (0 calls) | 0 (0 calls) |
| calls per data set: Null A / Null B | 5,775 / 1,466 | 2,573 / 0 | 7,358 / 1,287 | 5,926 / 644 |

What the table shows:
- **Null A** calls 27–75% false interactions on known-null data. Its hypothesis is false there, because
  one partner still depends on cell type, so this is a mismatch between the hypothesis and the claim,
  not a bug in Null A.
- **Null B** keeps FDR below 0.05 in all four scenarios. At α = 0.10 it stays below 0.10 (0.023–0.058).
- **Each half alone** fails exactly where its own side is untouched, so both halves are needed.
- **The family matters.** With BH over reported rows, Null B's FDR rises to 0.04–0.12. With BH over
  the testable rows, it is controlled.
- **Fully shuffled labels** (10 replicates per tissue) gave 0 calls by every procedure at every α,
  P(any call) ≤ 0.31 (Clopper–Pearson).

**Real labels (one run each, α = 0.05).**

| | Mammary | Liver |
|---|---|---|
| rows / testable / reported (prob > 0) | 518,144 / 71,322 / 8,952 | 1,251,648 / 80,850 / 15,839 |
| Null A calls | 6,876 | 12,753 |
| Null B calls | 3,290 (48% of Null A's) | 7,559 (59%) |
| called by Null B but not Null A | 0 | 0 |
| Null A only: the side that fails Null B (sender / receiver) | 1,297 / 2,289 | 3,216 / 1,978 |
| median prob, called by both / Null A only | 0.0071 / 0.0062 | 0.0089 / 0.0053 |
| Null B, δ = 0.01 / 0.05 | 0 / 0 | 1,781 / 0 |

- Null B keeps about half of Null A's calls. Every call it makes, Null A also makes.
- Each call it drops has one side that is no more specific to its cell type than a random labelling
  would make it.
- **The interval null is too strict for CellChat's scale.** prob is mostly below 0.01 (median of
  calls about 0.007), so δ = 0.01 is as large as the effect itself. **Frozen setting: δ = 0**, with
  BH over the testable rows.

**Resolution limit.** With B = 1,000 the smallest p is 1/1,001. At α = 0.01, BH over about 75,000
testable rows needs thousands of rows at that floor. Null B needs both sides at the floor, so at
α = 0.01 it made **0 calls in every scenario**, real labels included. Null A made 0 on Mammary and
12,287 on Liver. A B = 5,000 sensitivity run (Mammary, real labels and 3 receptor-shuffle replicates) is in progress and will be added here.

### Held-out test: TNBC (`null_b_tnbc.py`, run once with the frozen settings)

Open Problems `tnbc_wu2021`: 42,512 cells, 29 groups, 1,757,690 rows, 14,872 reported. B = 1,000,
seed 0, α = 0.05. `results/null_b_tnbc_wu2021.json`.

| | Null A (BH, testable) | Null B | Null B, δ = 0.01 |
|---|---|---|---|
| calls | 12,970 | 8,819 (68%) | 2,418 |
| Open Problems AUPRC (95% CI) | 0.350 (0.258–0.447) | 0.349 (0.254–0.459) | 0.303 (0.148–0.494) |
| Null B − this, AUPRC (95% CI) | −0.000 (−0.032 to 0.032) | | +0.046 (−0.101 to 0.189) |
| odds ratio | 3.45 | 3.45 | 2.20 |
| Jaccard of calls, patient halves (5 + 5 patients) | 0.38 (10,962 / 11,155 calls) | 0.30 (7,164 / 6,910) | 0.33 |

CellChat's unadjusted p < 0.05, and Null A over reported rows, give the same AUPRC (0.350) and
Jaccard (0.38) as Null A. Random scores give AUPRC 0.214. δ = 0.05 makes no calls, so its metrics
are not interpretable.

Calls of Null A only (4,151) against calls of both (8,819), on rows whose (ligand, target cell type)
the cytokine truth covers:

| | rows covered | truth-positive | precision (Clopper–Pearson 95%) |
|---|---|---|---|
| called by both | 153 | 57 | 0.37 (0.30–0.45) |
| Null A only | 32 | 15 | 0.47 (0.29–0.65) |

What TNBC shows:
- **No gain against the cytokine truth.** Null B removes a third of Null A's calls, and the AUPRC
  difference is 0.000 (CI ±0.03). This truth is "a ligand changes the target cell type's
  transcription", scored per (ligand, target) after aggregating over senders and receptors. It cannot
  tell whether the sender side was specific to s, which is half of what Null B tests. The rows Null B
  drops are not less often truth-positive (0.47 vs 0.37; overlapping CIs, 32 rows).
- **Lower reproducibility across patients** (Jaccard 0.30 vs 0.38). Part of this is the smaller call
  set: more of the calls sit near the threshold. Part is the max-p construction, since a call needs
  both sides significant in both halves. Null B's calls are not more stable than Null A's.
- Null B again made no call that Null A did not make.

**Conclusion.** Null B does what it was designed for. On data where "s signals to t" is known to be
false and one partner is still cell-type specific, it controls FDR (≤ 0.031), while Null A has 27–75%
false discoveries. On TNBC this stricter claim gives no measurable gain against the available truth,
and the call set is less reproducible across patients. Null B is the right test when the claim is a
specific sender-to-receiver edge. It does not improve recovery of ligand-to-target responses, and the
TNBC truth cannot judge the edge claim.

### What Null B does not answer

- **The receptor is not just "present" in t; it is specific to t.** A receptor expressed equally in
  every cell type fails H0_R, even if t uses it. Null B therefore tests the specificity of the
  s → t edge against other cell types, which is what CellChat's ranking is read as.
- **Candidate selection still uses the observed labels.** CellChat's over-expression filter
  (`identifyOverExpressedGenes`) picks genes by differential expression on the same labels. That step
  is the same for Null A and Null B, and its effect is the fully shuffled check above (0 calls).
- **Cells are the unit, not donors.** Permuting labels within donor (or permuting donors) would test
  reproducibility across individuals. That needs a multi-donor design and is future work.
- **Power at small α** is limited by B (see the resolution limit above).

## Checklist

- [x] cite one exact CellChat version: release 2.1.2 (plus main and the paper's version)
- [x] keep all rows (prob = 0 rows get p = 1), not only significant ones, so nulls can be compared
- [x] complexes keep CellChat's names in `ligand` / `receptor`; genes are in `ligand_complex` / `receptor_complex`
- [x] `nulls.null_a`: reproduce CellChat's own permutation test; compare with `p_native`
- [x] calibration: 20 shuffles, two tissues, CIs; ground truth with bootstrap CIs
- [x] `nulls.null_b`: composite null, calibrated on known-null data (20 replicates × 2 scenarios × 2 tissues), frozen (δ = 0, testable family), then run once on TNBC
- [x] B1 (empirical null) written up as a negative result

## Run the tests and the checks

```bash
pytest tools/cellchat
python -m tools.cellchat.match_null_a --quick
python -m tools.cellchat.calibrate_null_a --download --tissue all
python -m tools.cellchat.calibrate_null_a --tissue all
python -m tools.cellchat.null_b_calibration --tissue mammary --scenario receptor --reps 20   # also: ligand, labels, real; liver
python -m tools.cellchat.null_b_calibration --summarise
python -m tools.cellchat.null_b_tnbc
python -m tools.cellchat.b1_diagnosis --tissue all
python -m tools.cellchat.b1_tail_accuracy
```
