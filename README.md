# ccc-nulls

[![Tests](https://github.com/kshitij-bhardwaj/ccc-nulls/actions/workflows/tests.yml/badge.svg)](https://github.com/kshitij-bhardwaj/ccc-nulls/actions/workflows/tests.yml)
[![Docs](https://github.com/kshitij-bhardwaj/ccc-nulls/actions/workflows/docs.yml/badge.svg)](https://kshitij-bhardwaj.github.io/ccc-nulls/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.10%20%7C%203.11%20%7C%203.12-blue)
[![Kaggle runs](https://img.shields.io/badge/Kaggle-benchmark%20runs-20BEFF?logo=kaggle&logoColor=white)](kaggle/README.md)

**Is the sender specific *and* the receiver specific?**

Permutation-based cell–cell communication tools, such as CellPhoneDB, CellChat and the
[LIANA+](https://github.com/saezlab/liana-py) consensus built on them, test interactions under one null. They permute
every cell's label at once. That null is rejected as soon as **either** side is specific to its cell type. A ligand
made only by the sender, paired with a receptor every cell type expresses, is then called for every receiver.

`ccc-nulls` replaces their p-value with **Null B**, a composite null (intersection–union test). Under Null B the
ligand must be specific to the sender **and** the receptor must be specific to the receiver. It is computed from
LIANA+'s own permutations, so it needs no extra ones, and it slots into LIANA+'s consensus ranking. The output has
LIANA+'s columns, so LIANA+'s plots and downstream tools work unchanged.

## Installation

```bash
pip install "ccc-nulls[liana] @ git+https://github.com/kshitij-bhardwaj/ccc-nulls.git"
```

Python ≥ 3.10. The exact environment of the published results is in
[`requirements/lock-liana-final.txt`](requirements/lock-liana-final.txt) (LIANA+ 1.8.0).

## Quickstart

```python
import ccc_nulls as cn
import scanpy as sc

adata = sc.datasets.pbmc68k_reduced()
cn.rank_aggregate(adata, groupby="bulk_labels", use_raw=True)   # LIANA+'s rank_aggregate, with Null B
adata.uns["liana_res"].head()
```

| Argument | Default | |
|---|---|---|
| `null` | `"null_b"` | `"null_a"` gives LIANA+'s consensus with CellChat as a sixth method |
| `mode` | `"swap"` | `"swap"` replaces CellPhoneDB's and CellChat's p with Null B; `"extra"` adds Null B as two more methods |
| `strata` | `"auto"` | with several donors (`donor_id`, `patient`, `orig.ident`, … in `obs`), every permutation stays within a donor; `None`: global, as LIANA+ ([why](docs/method.md#donor-structure)) |
| `return_tools` | `False` | also return the CellPhoneDB and CellChat tables, with `p_L`, `p_R` and Null A and Null B p-values |

The [quickstart notebook](docs/notebooks/quickstart.ipynb) runs this end to end. On these PBMCs, 197 of the 437
rows that LIANA+'s CellPhoneDB calls at p < 0.05 are one-sided.

## How it works

For a candidate $L_s \to R_t$ with the tool's score $S$, permute the sender side only, then the receiver side only,
using LIANA+'s permutations:

$$p_L = \frac{1 + \#\{b : S(L_s^{\pi_b}, R_t) \ge S\}}{B+1},\qquad
p_R = \frac{1 + \#\{b : S(L_s, R_t^{\pi_b}) \ge S\}}{B+1},\qquad p_B = \max(p_L, p_R).$$

$p_B$ is a valid level-α test of "either side is unspecific", with no multiplicity correction (Berger 1982). The
magnitude scores and the rank aggregation are LIANA+'s own. Details: [Method](docs/method.md).

## Benchmark

Two pre-registered benchmarks answer two questions, in this order. Each has its own plan and amendments, every
one hashed before its runs.

**1. Are each tool's own p-values valid?** ([`benchmarks/liana_nulls`](benchmarks/liana_nulls/PLAN.md);
results: [`benchmarks/results/liana_nulls`](benchmarks/results/liana_nulls/README.md))
- **Three columns per tool:** the tool itself (CellChat 2.1.2 in R, CellPhoneDB 5.0.1), our Python Null A on the
  tool's own permutations, and Null B. Columns 1 and 2 agree on every row of five data sets.
- **Calibration:** on Tabula Sapiens with one side's genes shuffled, Null A's FDR at 5% is 0.27–1.00; Null B's
  is 0–0.03.
- **Power:** with BH, Null B needs about 5,000 permutations (CellChat at f = 1.5: power 0.19 at 1,000, 0.97 at
  5,000).
- **Also compared:** B2 (decoy genes, CAMERA-corrected) and PR #9's CellPhoneDB Null B.

**2. Does Null B improve LIANA+'s consensus ranking?** ([`benchmarks/liana_final`](benchmarks/liana_final/PLAN.md),
then [`AMENDMENT_1.md`](benchmarks/liana_final/AMENDMENT_1.md) to
[`AMENDMENT_4.md`](benchmarks/liana_final/AMENDMENT_4.md); results:
[`benchmarks/results/liana_final`](benchmarks/results/liana_final/README.md))
- **What is compared:** LIANA+'s six methods (its five plus CellChat), with the p-values of CellPhoneDB and CellChat
  replaced by Null B ("swap") or added as two more methods ("extra"), against LIANA+ as it is.
- **Measures:** calibration, power and biology (spatial and downstream-response truths), then a held-out
  confirmation on Tabula Sapiens Lung, donor-aware permutations, sensitivity, a second species and comparison with
  simple filters (AMENDMENT_2).
- **Second candidate:** AMENDMENTS 3 and 4 add T2 (PR #11: a leave-one-out test of each cell-type pair against the
  interaction's other pairs) to every run.

Both cover the two Open Problems communication sets that LIANA+'s authors built (TNBC and Allen brain), Tabula
Sapiens Mammary and Liver, and the mouse hypothalamus; the second adds Tabula Sapiens Lung. Full results:
[Benchmark](docs/benchmark.md).

The second benchmark runs on Kaggle (4 CPUs, 32 GB per session), one session per part, from one commit and the locked
environment. Each result file records the link of the session that produced it, and
[`kaggle/README.md`](kaggle/README.md) lists every session with its link, inputs, command and run time. The inputs
are three Kaggle data sets: the [code](https://www.kaggle.com/datasets/kshitijbhardwaj1729/ccc-nulls-code) (a `git archive` of the commit), the [auxiliary files](https://www.kaggle.com/datasets/kshitijbhardwaj1729/ccc-nulls-aux)
(ligand–receptor databases, CytoSig, PROGENy, HCOP) and the [Open Problems copies](https://www.kaggle.com/datasets/kshitijbhardwaj1729/ccc-nulls-openproblems). The sessions and data sets
are private to the project team; to be added as a collaborator, ask a maintainer.

<!-- RESULTS -->

## Reproducing

```bash
make install    # the locked environment
make data       # ~9 GB, md5-checked
make reproduce  # every result in benchmarks/results/liana_final/
```

To rerun everything on Kaggle, see [`kaggle/README.md`](kaggle/README.md). To split the runs across machines or run over SSH, see [`RUN.md`](RUN.md). Every part resumes where it stopped. Every result file records the SHA-256 of the plan and amendments it was
produced under. See [Reproducibility](docs/reproducibility.md).

## Repository layout

```
src/ccc_nulls/          the package: rank_aggregate, liana_consensus, composite, competitive (B2), datasets
benchmarks/liana_nulls/ benchmark 1, each tool's own p-values: PLAN, amendments, code
benchmarks/liana_final/ benchmark 2, LIANA+'s consensus: PLAN, amendments, code
benchmarks/results/     result JSON files, results READMEs, figures
docs/                   documentation (Sphinx) and the quickstart notebook
tools/                  native adapters for CellPhoneDB, CellChat, NATMI, Connectome
tests/                  tests
```

The earlier tool-agnostic significance layer (`ccc_nulls.significance_layer`) and its simulation and TNBC results
are documented in [Earlier work](docs/earlier_work.md).

## Citation

If you use ccc-nulls, cite it with [`CITATION.cff`](CITATION.cff), and cite LIANA+, which it extends:

> Dimitrov D. et al. LIANA+ provides an all-in-one framework for cell–cell communication inference.
> *Nat Cell Biol* 26, 1613–1622 (2024).

## License

[BSD-3-Clause](LICENSE), as LIANA+'s.
