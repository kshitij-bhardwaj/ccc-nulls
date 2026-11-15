# Null A vs Null B inside the LIANA+ pipeline: analysis plan

Fixed on 1 November 2026, before any result below was computed. Every result file records the
SHA-256 of this file. If the plan changes, the change and its reason are added at the end, with a
new date. Earlier numbers are never edited.

## Question

Each tool's own permutation null (Null A) asks whether *either* side of a ligand–receptor edge is
specific to its cell type. A call is read as "cell type s signals to cell type t through L–R".
Null B (`ccc_nulls.composite`) tests that claim directly: the ligand side must be specific to s
**and** the receptor side to t (an intersection–union test, p = max(p_L, p_R)).

We ask:
- **validity:** which null controls false discoveries for the s → t claim?
- **power:** what does Null B lose?
- **biology:** is the edge set Null B keeps more biologically plausible than the one Null A keeps?

## Pipeline (one script, `python -m benchmarks.liana_nulls`)

For each dataset and tool (CellChat 2.1.2, CellPhoneDB 5.0.1), through the tool's adapter:

1. **LIANA+ as it is.** `li.mt.rank_aggregate` plus the tool's own LIANA+ method; their output is
   kept unchanged.
2. **Column 1, native.** The tool itself with its defaults:
   - CellChat in R: `computeCommunProb`, nboot = 100, seed.use = 1;
   - CellPhoneDB statistical analysis: 1,000 iterations, threshold 0.1, debug_seed = 0, threads = 1.
   - The tool's call rule: p < 0.05 (CellChat), p ≤ 0.05 (CellPhoneDB `significant_means`).
3. **Column 2, our Python Null A,** on the same permutation sequence as the tool. Acceptance
   criterion: score, p-value and call set identical to column 1 on every row.
4. **Column 3, Null B,** through the adapter's `side_score`:
   - the tool's exact score, the same permutations as Null A, B = 1,000;
   - BH at α = 0.05 over the same rows as Null A (the family frozen in each tool's calibration).
5. **Extra column:** LIANA+'s own re-implementation of the tool. It is not expected to equal
   column 1; the agreement is reported.

Once column 1 and column 2 are shown identical on a dataset, every later Null A number on that
dataset is computed with column 2. It is then also the native tool's number.

## Stages and endpoints

| Stage | Data | Endpoint (fixed) |
|---|---|---|
| T0 fidelity | TNBC, Tabula Sapiens Mammary, Liver | # rows where column 1 ≠ column 2 (must be 0) |
| T1 calibration | Tabula Sapiens Mammary, Liver | FDR = E[V/R] on known-null data (one side shuffled), mean with t-CI over 20 replicates; P(any call) on shuffled labels |
| T2 power | Tabula Sapiens Mammary, labels shuffled, edges planted (below) | FDR and power on planted edges, per effect size, mean with t-CI over 5 replicates |
| T3 biology | TNBC (Wu 2021 Visium truth; Open Problems cytokine truth), Allen brain (Open Problems spatial truth) | see below |

**T2 planted edges.**
- Background: Tabula Sapiens Mammary, 10,000 cells, cell-type labels shuffled, so no edge is true.
- Each replicate plants 30 database L–R pairs, each with a random (s, t), s ≠ t:
  - the ligand is raised in s and the receptor in t;
  - genes are untouched elsewhere.
- Raising a gene by fold change f in a cell type:
  - each non-zero value becomes log1p(f · expm1(x));
  - each zero is switched on with probability 1 − 1/f, taking a value drawn from that gene's non-zero values in the type, raised the same way.
- Effect sizes f = 1.5, 2, 3, 4.
- True rows: the planted (s, t, L, R). Every other row is null for the s → t claim, including:
  - (s, t′, L, R) with t′ ≠ t: ligand specific, receptor not specific to t′;
  - (s′, t, L, R) with s′ ≠ s.
- Reported per method: power = true calls / planted rows, FDP, and the false calls that sit on a planted pair (one-sided).

**T3 truths and metrics.** Each truth is labelled by which side of the edge it can see:
- **Spatial colocalisation (sees s and t).**
  - TNBC: the 4 TNBC Visium slides of Wu et al. 2021 (Zenodo 4739739).
    - Spot deconvolution onto the 29 atlas cell types by non-negative least squares on cell-type mean profiles of marker genes (SPOTlight's final step).
    - The method follows Dimitrov et al. 2022: per slide, the Pearson correlation of cell-type proportions across spots; z-scores over all s ≠ t pairs averaged over slides; z ≥ 1.645 means colocalised.
  - Allen brain: Open Problems' own source–target truth, used as given.
- **Cytokine activity (sees L and t):** the Open Problems truth on TNBC, as before.
- **Fake gene pairs (sees whether L and R interact):** as before.

Metrics:
- **Primary (spatial, per dataset):** among Null A's calls with s ≠ t, compare the calls Null B keeps with the calls it drops.
  - Colocalisation odds ratio, kept vs dropped.
  - 95% CI by bootstrap over (s, t) pairs, 2,000 resamples, because calls share (s, t).
- **Secondary:**
  - (s, t)-level AUPRC, scoring each pair by its number of calls, with a paired bootstrap difference Null B − Null A;
  - the cytokine AUPRC difference;
  - the fake-pair share of calls.

## Decision rule (fixed)

- **Null B is the better null** for the s → t claim if all three hold:
  1. its T1 FDR upper CI bound is ≤ 0.05 in every known-null scenario, and its T2 FDR upper bound is ≤ 0.05 at every effect size;
  2. the primary spatial odds ratio (kept vs dropped) has a 95% CI above 1 in at least one spatial dataset and below 1 in none;
  3. no secondary truth shows Null B significantly worse: a paired-bootstrap CI for its AUPRC entirely below −0.02.
- **Valid but not biologically distinguishable:** (1) and (3) hold, but (2) does not.
- **Rejected:** (1) fails, or (3) fails.

Null A is judged by the same criteria, so the two can be compared directly.

## Not in this plan

Patient- and donor-stratified permutations, and patient-split reproducibility with CIs, come later.
