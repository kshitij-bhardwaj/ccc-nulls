# ccc-nulls: code

Permutation-based cell–cell communication tools (CellPhoneDB, CellChat, and LIANA+'s consensus built on them)
call an interaction when **either** the ligand side or the receptor side is specific to its cell type. This
package adds tests whose null is "not both sides specific": the bilateral permutation test (Null B in the code,
p = max(p_L, p_R) from one-sided label permutations), a matched-decoy test and a leave-one-out t-test. Two
pre-registered benchmarks compare them with each tool's own test: at tool level (Benchmark 1) and inside LIANA+'s
consensus ranking (Benchmark 2). Every number in the report comes from the JSON files in `benchmarks/results/` and
`tools/cellchat/results/`, produced by the commands below.

Names in the report and slides vs. names in the code and result files:

| Report / slides | Code and result files |
|---|---|
| Joint label permutation (JLP), the tools' own test | Null A, `null_a` |
| Bilateral permutation test (BPT); stratified BPT = within-donor | Null B, `null_b` (`ccc_nulls.composite`, `ccc_nulls.liana_consensus`) |
| Matched-decoy test (MDT) | B2, Null B2 (`ccc_nulls.competitive`) |
| Leave-one-out t-test (LOO-t) | T2 (`ccc_nulls.type_perm`, `ccc_nulls.nullb_type_restd`) |
| -R (replace the tool's p) / -A (add as extra methods) | `swap_*` / `extra_*` (7th/8th method) |
| LIANA+ baseline | `liana6` (LIANA+'s five methods plus CellChat) |

## 1. Installation

Python 3.12 (the package itself supports 3.10–3.12; the results were produced with 3.12). From this `code/` folder:

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt          # the locked environment of the results (LIANA+ 1.8.0, CellPhoneDB 5.0.1)
pip install --no-deps -e .               # the ccc_nulls package from src/
export PYTHONPATH=$PWD/src:$PWD          # the benchmarks and tools run as modules from this folder
```

(`make install` does the two `pip` steps.) Quick check, about 15 seconds, no download (scanpy ships these PBMCs):

```bash
python -c "import ccc_nulls as cn, scanpy as sc; a = sc.datasets.pbmc68k_reduced(); \
cn.rank_aggregate(a, groupby='bulk_labels', use_raw=True); print(a.uns['liana_res'].head())"
```

`cn.rank_aggregate` is LIANA+'s `rank_aggregate` with CellPhoneDB's and CellChat's p-values replaced by Null B;
its output has LIANA+'s columns.

**R (only for Benchmark 1's native CellChat column).** R ≥ 4.3 with CellChat 2.1.2:

```bash
Rscript tools/cellchat/install_r_packages.R      # or: make install-r
```

The R session of the native CellChat runs (R 4.3.2, CellChat 2.1.2, Matrix 1.6-1.1 and every loaded package) is in
[`tools/cellchat/R_sessionInfo.txt`](tools/cellchat/R_sessionInfo.txt). Nothing else uses R: the null tests run on
our Python port of CellChat, which matches CellChat in R row for row
([`tools/cellchat/README.md`](tools/cellchat/README.md)).

**Native NATMI (optional, Benchmark 1 AMENDMENT_5's column 1 only)** runs in its own Python 3.7 environment;
[`tools/natmi/requirements-native.txt`](tools/natmi/requirements-native.txt) has the steps.

**Tests** (no downloaded data; the CellChat tests fetch CellChatDB from GitHub once):

```bash
python -m pytest -q        # or: make test
```

## 2. Data

No data is in the repository. Data go to `data/` (or the folder in `$CCC_NULLS_DATA`); file names and md5 sums are
in [`data/README.md`](data/README.md).

```bash
python -m ccc_nulls.datasets tnbc_wu2021 allen_brain_atlas tabula_sapiens_mammary tabula_sapiens_liver tabula_sapiens_lung
```

This downloads in 16 parallel byte ranges and checks each file's md5 (about 9 GB in all; `make data`).

| Data set | Role | Source |
|---|---|---|
| Tabula Sapiens Mammary, Liver (human) | calibration and power (known-null shuffles, planted pairs) | CELLxGENE collection https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5; files https://datasets.cellxgene.cziscience.com/33b25f68-a9d6-4643-9283-653e937ae608.h5ad (Mammary), https://datasets.cellxgene.cziscience.com/a4b3a49e-062b-4e3f-8915-02f40607651f.h5ad (Liver) |
| Tabula Sapiens Lung (human) | held-out confirmation (Benchmark 2, AMENDMENT_2 on) | https://datasets.cellxgene.cziscience.com/fddf43a3-b7c8-49e9-ac0a-5e24b390849f.h5ad |
| TNBC, Wu et al. 2021 (human, 10 patients) | real data; cytokine truth; donor tests | Open Problems: https://openproblems-data.s3.amazonaws.com/resources/datasets/openproblems_v1/tnbc_wu2021/log_cp10k/dataset.h5ad |
| Allen brain, mouse V1 | real data; spatial truth (Open Problems' own) | Open Problems: https://openproblems-data.s3.amazonaws.com/resources/datasets/openproblems_v1/allen_brain_atlas/log_cp10k/dataset.h5ad |
| Mouse hypothalamus, Moffitt et al. 2018 | real data; MERFISH spatial truth | scRNA-seq: GEO GSE113576 (the 10x files `GSE113576_{matrix.mtx,genes.tsv,barcodes.tsv}.gz` in `data/moffitt2018/scrna/`); MERFISH: animal 1 of the study (scverse example data) as `data/moffitt2018/merfish/merfish.h5ad`. `benchmarks/liana_nulls/hypothalamus.py` builds `data/moffitt2018/moffitt2018_hypothalamus.h5ad` from them |
| TNBC Visium slides, Wu et al. 2021 | TNBC spatial truth | Zenodo record 4739739, the 4 TNBC slides, in `data/wu2021_visium/` |

**Open Problems closed its S3 bucket** (`AllAccessDisabled`, checked 2026-11-05), so the two Open Problems links no
longer download. Copies taken before that are the Kaggle data set
[`kshitijbhardwaj1729/ccc-nulls-openproblems`](https://www.kaggle.com/datasets/kshitijbhardwaj1729/ccc-nulls-openproblems)
(private to the project team: ask to be added as a collaborator). Put the two files in `data/` under the names in
`data/README.md`; the command above then checks their md5 against the originals.

**Auxiliary files** (about 200 MB): the CellChatDB v2.1.2 files (also fetched automatically from
`github.com/jinworks/CellChat` the first time), the CellPhoneDB v5.0.0 database (`cellphonedb_v5.0.0/cellphonedb.zip`),
the CytoSig signature (`downstream/cytosig_signature.centroid`, Jiang et al. 2021), PROGENy from OmniPath
(`downstream/progeny_omnipath_human.tsv`, Schubert et al. 2018), HCOP human–mouse orthologs (also fetched through
LIANA+ the first time) and the processed hypothalamus. They are one archive, the Kaggle data set
[`kshitijbhardwaj1729/ccc-nulls-aux`](https://www.kaggle.com/datasets/kshitijbhardwaj1729/ccc-nulls-aux) (private, as above):

```bash
tar xzf ccc-nulls-aux.tgz -C data
(cd data && md5sum -c ../benchmarks/liana_final/aux_files.md5)      # macOS: brew install coreutils
```

The spatial truths derived from the Visium and MERFISH data are committed
(`benchmarks/results/liana_nulls/truth_*_spatial.{csv,json}`), so the raw spatial files are needed only to recompute
them (`python -m benchmarks.liana_nulls.spatial`).

## 3. Reproducing the results

Every part writes its JSON after each run and skips runs already written, so a stopped part resumes where it
stopped. Every result file records the SHA-256 of the plan and amendment it was produced under.

### Results tables, from the committed JSON files (seconds, no data)

```bash
make reports
```

This regenerates every results README (the tables quoted in the report) and the Benchmark 2 summary figures in
`benchmarks/results/liana_final/figures/`. The commands it runs:

| Output | Command |
|---|---|
| `benchmarks/results/liana_nulls/README.md` (Benchmark 1, all amendments) | `python -m benchmarks.liana_nulls --part report` |
| `benchmarks/results/liana_final/README.md` (Benchmark 2, PLAN + AMENDMENT_1) | `python -m benchmarks.liana_final --part report` |
| `benchmarks/results/liana_final/amendment<N>/README.md`, N = 2, 3, 4, 6, 8 | `python -m benchmarks.liana_final.amendment<N> --part report` |
| `benchmarks/results/liana_final/amendment7/README.md` | `python -m benchmarks.liana_final.report_amendment7` |
| `benchmarks/results/liana_final/figures/*.pdf` | `python -m benchmarks.liana_final.figures` |

AMENDMENT_5 of Benchmark 2 has no README of its own: the AMENDMENT_6 and AMENDMENT_7 reports read its JSON files.

### Benchmark 1: is each tool's own p-value valid? (`benchmarks/liana_nulls`)

Plan: [`benchmarks/liana_nulls/PLAN.md`](benchmarks/liana_nulls/PLAN.md), then `AMENDMENT_1.md` … `AMENDMENT_5.md`.
Three columns per tool: the tool itself (CellChat 2.1.2 in R, CellPhoneDB 5.0.1), our Python Null A on the tool's
own permutations, and Null B; plus B2 and T2. Needs the data, the auxiliary files and R. Results:
`benchmarks/results/liana_nulls/`.

```bash
python -m benchmarks.liana_nulls --part core                                      # T0: columns 1-3 on TNBC, Allen brain, Mammary, Liver
python -m benchmarks.liana_nulls --part core --dataset moffitt2018_hypothalamus   # AMENDMENT_1's data set
python -m benchmarks.liana_nulls --part calibration power biology amendment extend power_a3 b2_validity ablation report
python -m benchmarks.liana_nulls.amendment4 --part all                            # T2 (LOO-t) through every stage
python -m benchmarks.liana_nulls.amendment5 --part all                            # NATMI (native NATMI for column 1)
```

(`make benchmark1` runs the first four.) The calibration stage reads each tool's calibration runs, or runs them
when their files are missing: `python -m tools.cellchat.null_b_calibration --tissue mammary` (and `liver`;
results in `tools/cellchat/results/`) and
`python -m benchmarks.null_b_cellphonedb --db data/cellphonedb_v5.0.0/cellphonedb.zip --part calibration --tissue mammary`
(and `liver`; results in `benchmarks/results/null_b_cellphonedb/`). The CellChat port's own checks against CellChat
are `python -m tools.cellchat.match_null_a` and `python -m tools.cellchat.calibrate_null_a --tissue all`.

### Benchmark 2: does Null B improve LIANA+'s consensus? (`benchmarks/liana_final`)

Plan: [`benchmarks/liana_final/PLAN.md`](benchmarks/liana_final/PLAN.md), then `AMENDMENT_1.md` … `AMENDMENT_8.md`:
held-out Lung, donors, sensitivity, mouse and filters (2), T2 (3, 4), the same tests for every null (5), Connectome's
p replaced too (6, 7), NATMI's too (8). Results: `benchmarks/results/liana_final/` and its `amendment<N>/` folders.

```bash
python -m benchmarks.liana_final --part all              # PLAN + AMENDMENT_1: run, calibration, power, biology
python -m benchmarks.liana_final.amendment2 --part all
python -m benchmarks.liana_final.amendment3 --part all
python -m benchmarks.liana_final.amendment4 --part all
python -m benchmarks.liana_final.amendment5 --part all
python -m benchmarks.liana_final.amendment6 --part all
python -m benchmarks.liana_final.amendment7 --part all
python -m benchmarks.liana_final.amendment8 --part all
```

(`make benchmark2`.) This is days of computing on one machine with 32 GB of RAM. The published results ran as
Kaggle sessions (4 CPUs, 32 GB each), one per part, from one commit and `requirements.txt`; the runs that needed
more memory ran on a lab server. [`kaggle/README.md`](kaggle/README.md) lists every session with its exact command,
inputs and run time, and `kaggle/kernels.py` builds, dispatches and collects them
(`python kaggle/kernels.py code && python kaggle/kernels.py build && python kaggle/kernels.py dispatch && python kaggle/kernels.py collect`).
Each result file records the session or machine that produced it. To split the parts across your own machines,
see [`RUN.md`](RUN.md); `--only REGEX` keeps only the runs whose tag matches.

## 4. Paper figures

The report's three figures come from the committed result files; nothing is recomputed (a few seconds):

```bash
cd code
python figures/make_figures.py            # all three; or name some: python figures/make_figures.py fig2
```

| Output (`figures/output/`) | Report | What it shows |
|---|---|---|
| `fig1_tool_level.pdf` | Fig. 1 | Benchmark 1: FDR on known nulls and power vs FDR, all four tools |
| `fig2_consensus.pdf` | Fig. 2 (main figure) | Benchmark 2: LIANA+'s top 1% and planted-pair power as each test replaces 2, 3, 4 tools' p-values |
| `fig3_robustness.pdf` | Fig. 3 | p-values inside LIANA+, donor shuffles, sensitivity runs, spatial truth |

`figures/final_results.json` holds the PLAN's Benchmark-1 numbers that the figures read, copied from
`benchmarks/results/liana_nulls/README.md` and `benchmarks/results/liana_final/README.md` (its `source` field).

## 5. Folder structure

```
code/
├── README.md               this file
├── requirements.txt        the locked Python environment of the results
├── pyproject.toml          the ccc_nulls package (pip install -e .) and pytest settings
├── Makefile                install, install-r, test, data, benchmark1, benchmark2, reports
├── figures/                make_figures.py: the report's figures from the committed results
├── RUN.md                  splitting Benchmark 2 across machines, and merging the results
├── LICENSE                 BSD-3-Clause, as LIANA+
├── src/ccc_nulls/          the package
│   ├── liana_consensus.py    rank_aggregate: LIANA+'s consensus with Null B (or B2, T2) inside; within-donor by default
│   ├── composite.py          Null B: one-sided permutations, p = max(p_L, p_R)
│   ├── competitive.py        B2: matched decoy genes, CAMERA-corrected
│   ├── type_perm.py, nullb_type_restd.py                  T2: the leave-one-out t-test
│   ├── nulls.py, core.py, scores.py, stats.py, _matrix.py   Null A and the shared permutation and BH code
│   ├── liana_reference.py, schema.py, adapters.py, io.py    LIANA+'s tool outputs in one common schema
│   ├── datasets.py           data set registry, downloader with md5 check, loaders
│   ├── evaluation/           Open Problems' evaluation metrics
│   └── simulate.py           synthetic data with planted interactions (used by the tests)
├── benchmarks/
│   ├── liana_nulls/          Benchmark 1: PLAN.md, AMENDMENT_1-5.md and their code
│   ├── liana_final/          Benchmark 2: PLAN.md, AMENDMENT_1-8.md and their code (amendment<N>.py runs
│   │                         AMENDMENT_<N>.md, report_amendment<N>.py writes its README), aux_files.md5
│   ├── null_b_cellphonedb.py CellPhoneDB's own Null B (Benchmark 1's "#9's Null B" column, AMENDMENT_3)
│   └── results/              every result: JSON files, results READMEs (the tables), figures
├── tools/                  one folder per tool: native run, conversion to the common schema, Null A port
│   ├── cellchat/             CellChat 2.1.2 ported to Python; native R run (run_cellchat.R, native_r.py);
│   │                         R install script and sessionInfo; port checks and calibration (results/)
│   ├── cellphonedb/          CellPhoneDB 5.0.1 adapter, Null A and Null B on CellPhoneDB's statistic
│   └── natmi/                NATMI adapter (its native environment: requirements-native.txt)
├── kaggle/                 the Kaggle sessions that produced Benchmark 2 (kernels.py, runner.py, watch.py, README.md)
├── tests/                  tests of the package and the benchmark code
└── data/                   downloaded data (git-ignored); data/README.md lists every file and its md5
```

The plan documents (`PLAN.md`, `AMENDMENT_*.md`) were committed before the runs they govern, and every result file
records their SHA-256, so they are left exactly as written. They call the locked environment
`requirements/lock-liana-final.txt` (now `requirements.txt`) and refer to pull requests of the development
repository (#9: CellPhoneDB's Null B; #11: T2).
