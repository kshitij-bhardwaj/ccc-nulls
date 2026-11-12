# Amendment 5: NATMI as a third tool, through every stage and every null

Written on 2026-11-10, after every result of `PLAN.md` and Amendments 1–4 was known, and before any null had been
computed on NATMI on any data set. This is an **addition**: it changes no earlier number, verdict or result file.
Result files written under it are `amendment5_*.json`; each records the SHA-256 of `PLAN.md` and of this file.

**What had been seen when this was written:**
- every result of this benchmark (`benchmarks/results/liana_nulls/README.md`) and of `benchmarks/liana_final`;
- NATMI's source code (`asrhou/NATMI`, `ExtractEdges.py`, commit `f35f677cb2e4c8a176a6c5501a9483c9e28661a1`, the
  latest, 2021-04-13), to define its score;
- no NATMI output, native or ours, on any data set of this benchmark.

**Why.** The project brief asks for the significance layer on several tools. This benchmark ran two
(CellChat, CellPhoneDB). NATMI is the third: a widely used tool whose score has no permutation null at all.

Code: `tools/natmi/adapter.py` (the tool), `benchmarks/liana_nulls/tools.py` (`NATMI`),
`benchmarks/liana_nulls/amendment5.py` (every stage). Run: `python -m benchmarks.liana_nulls.amendment5`.

## The tool

- **Native NATMI:** `ExtractEdges.py` at the commit above, `--interDB lrc2p` (connectomeDB2020, literature-supported
  pairs), `--species human --interSpecies human --idType symbol`, in its own environment
  (`tools/natmi/requirements-native.txt`: the versions NATMI was tested with).
- **Score:** NATMI's *edge average expression derived specificity* (LIANA+'s `spec_weight`):

  spec(L_s → R_t) = (mean of L in s / Σ_k mean of L in k) · (mean of R in t / Σ_k mean of R in k),

  the means over all cells of the type (zeros included), Σ_k over every cell type; 0 where a sum is 0. NATMI's
  edge average expression weight (`expr_prod`, mean · mean) is recorded on every row but is not tested.
- **Input:** the benchmark's X = ln(1 + CP10k), as for the other two tools (NATMI's documentation suggests a linear
  scale; the benchmark's scale keeps the three tools on the same data). Cell types as for the other tools.
- **Mouse data** (Allen brain, hypothalamus): genes renamed to their HCOP one-to-one human orthologs first, as
  for CellPhoneDB (`data.to_human`), then run as human.
- **Rows:** every lrc2p pair with both genes in the data × every (source, target), including score-0 rows, built by
  `tools/natmi/adapter.py`. NATMI itself writes only rows with a positive score; those are the rows the fidelity
  check compares.

## The three columns

NATMI has no p-value, so its columns differ from the other two tools in one way:

| Column | CellChat, CellPhoneDB | NATMI |
|---|---|---|
| 1, native | the tool's own p | NATMI's own **score** (no p exists) |
| 2, our Python | Null A on the tool's own permutations | NATMI's score recomputed in Python, and Null A on it |
| 3 | Null B | Null B |

- **T0 fidelity (NATMI):** native and Python rows with a positive score must be the same set, and their scores
  must agree to a relative difference of 1e-9 (NATMI writes and re-reads its edges through `.xlsx` files, so exact
  float equality is not expected). If T0 fails on a data set, no NATMI number is reported for it.
- **Permutations:** NATMI has no permutation sequence of its own. Null A and Null B use the same
  `numpy.random.default_rng(seed)` sequence (one `permutation` of the cell-type codes per iteration), seed 0,
  B = 1,000, ties counted (p = (1 + #{s_b ≥ s}) / (B + 1)), the repository's convention.
- **Family (BH):** `composite.testable`: rows where some labelling (observed, Null A's or either side of Null B's)
  gives a positive score, as for CellChat.
- **The tool's own rule (Amendment 2 analogue):** NATMI has none. Its stand-in is CellChat's form: Null A
  p < 0.05 and score > 0, unadjusted. It is reported as `native_rule` / `tool_rule` and is not a NATMI output.
- **LIANA+'s re-implementation** (the extra column): `li.mt.natmi` with `resource_name="connectomedb2020"`.

## The nulls

| Null | On NATMI |
|---|---|
| Null A | as above |
| Null B | `ccc_nulls.composite` on NATMI's score (`tools/natmi/adapter.side_score`): the sender side is the ligand's specificity, the receiver side the receptor's; p = max(p_L, p_R) |
| B2 (MDT) | `ccc_nulls.competitive`, unchanged, on NATMI's rows and family. Decoys exclude every gene of every database already excluded **and** every lrc2p gene. B2 does not read the tool's score, so on rows NATMI shares with another tool its p equals that tool's B2 p |
| LOO-t | `ccc_nulls.nullb_type_restd.p_values` on NATMI's score over NATMI's family, as Amendment 4 |

## Stages (data and seeds as for the other two tools)

| Stage | As in | Runs | Reported |
|---|---|---|---|
| core (T0) | the plan's T0, Amendment 3's columns | TNBC, Allen brain, Tabula Sapiens Mammary and Liver, hypothalamus | fidelity; calls of every null (BH 0.05); Null B outside Null A; Null A only by limiting partner; LIANA+'s NATMI vs ours |
| calibration | `tools/cellchat/null_b_calibration.py`'s design | Tabula Sapiens Mammary and Liver; 20 replicates with receptor-only genes shuffled (seeds 2000 + rep), 20 with ligand-only genes shuffled (seeds 3000 + rep), 10 with labels shuffled (seeds 1000 + rep); lrc2p's ligand-only and receptor-only genes | FDR at BH 0.05 (mean, 95% t-CI) on rows whose shuffled side is null; P(any call) with labels shuffled; for Null A, Null B, B2 and LOO-t on the same replicates |
| power (T2) | the plan's T2 | the plan's planted data unchanged (replicates 0–4, f = 1.5, 2, 3, 4) | power and FDR (mean, 95% t-CI) of every null, BH 0.05, plus the tool-rule stand-in; power is over the planted pairs in lrc2p, and their number is reported |
| biology (T3) | the plan's T3 and Amendment 1 | spatial: TNBC Visium, Allen brain, hypothalamus; cytokine (TNBC); CytoSig and PROGENy on all five | odds ratio of colocalisation (Null A's calls kept vs dropped, 2,000 bootstrap resamples of (s, t)) for Null B, B2 and LOO-t; pair-level AUPRC; D1 for the downstream truths |

The calibration harness differs from the one the plan used for the other two tools (separate scripts per tool,
and Amendments 3 and 4's 10 replicates with seeds 4000–6000 for B2 and LOO-t): here all four nulls are scored on
the same replicates. The false rows are defined the same way: rows whose shuffled side has every gene shuffled.
As in CellChat's calibration script, replicate r's permutations use seed r + 1 (Null A and Null B alike), and,
as in Amendments 3 and 4, B2 is tested on the core rows and family of the unshuffled tissue.

Tabula Sapiens is loaded with lrc2p's genes added to the genes the benchmark already loads (and B2's decoy
universe); library sizes use every gene, so no value of the other tools' genes changes. The power stage's
candidate pairs, backgrounds and planted genes are computed from the same genes as before and are unchanged.

## Decision

The plan's rule, applied to NATMI's Null A and Null B as written, and to B2 and LOO-t as Amendments 3 and 4 apply
it: **valid** if the FDR's upper 95% bound is ≤ 0.05 in every calibration setting and in the power stage; the
biology criteria as in the plan. NATMI's results are reported next to the other two tools, never pooled with them.
