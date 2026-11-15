# CS690A Project 6: structure-preserving null models and FDR control for cell–cell communication

Permutation-based cell–cell communication tools (CellPhoneDB, CellChat, and LIANA+'s consensus built on them)
test each ligand–receptor interaction by permuting every cell's label at once, so an interaction is called as
soon as either the ligand or the receptor is specific to its cell type. We test instead whether both sides are
specific, with a bilateral permutation test (BPT, computed from the tools' own permutations), and compare it with a
matched-decoy test, a leave-one-out t-test and the tools' own test in two pre-registered benchmarks: each tool's
p-values, and LIANA+'s consensus ranking, on Tabula Sapiens, TNBC, Allen brain and mouse hypothalamus data with
known-null shuffles, planted interactions and spatial and downstream truths.

## Layout

```
CS690_Project6/
├── README.md     this file
├── slides/       final presentation: PDF and source
├── report/       report.pdf
└── code/         the final code, with its own README.md: installation (Python and R), data sources and download,
                  the command for each result, folder structure
```

To reproduce: follow [`code/README.md`](code/README.md); its section 5 names the code behind every table and figure. In short, from `code/`:

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt && pip install --no-deps -e .
export PYTHONPATH=$PWD/src:$PWD
python -m pytest -q          # tests, no data needed
make reports                 # every results table from the committed result files (seconds)
python figures/make_figures.py   # the report's three figures (Fig. 2 is the main one) from the same files
make data                    # the data sets (about 9 GB), then: make benchmark1, make benchmark2
```
