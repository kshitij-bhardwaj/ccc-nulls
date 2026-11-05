# Running the benchmark

Short guide for running `benchmarks/liana_final` on your own machine or over SSH. Every part resumes where it
stopped, so a killed run loses at most the run in progress.

## 1. Set up (Python 3.12, about 5 minutes)

```bash
git clone -b benchmarks/publication https://github.com/kshitij-bhardwaj/ccc-nulls.git && cd ccc-nulls
uv venv --python 3.12 .venv && source .venv/bin/activate        # or: python3.12 -m venv .venv
uv pip install -r requirements/lock-liana-final.txt && uv pip install --no-deps -e .   # plain pip works too
export CCC_NULLS_DATA=$PWD/data PYTHONPATH=$PWD                    # put both in ~/.bashrc on a server
```

## 2. Get the data (about 10 GB)

```bash
python -m ccc_nulls.datasets tnbc_wu2021 allen_brain_atlas tabula_sapiens_mammary tabula_sapiens_liver tabula_sapiens_lung
```

The auxiliary inputs (CellChat and CellPhoneDB databases, CytoSig and PROGENy, HCOP, the processed
hypothalamus) are one archive, `ccc-nulls-aux.tgz` (200 MB). It is the Kaggle dataset `kshitijbhardwaj1729/ccc-nulls-aux`,
which is private: ask to be added, or copy it from a teammate. Then:

```bash
tar xzf ccc-nulls-aux.tgz -C data
(cd data && md5sum -c ../benchmarks/liana_final/aux_files.md5)   # macOS: brew install coreutils
```

## 3. Run

Everything, on one machine with at least 32 GB of RAM (about 1.5 days on 8 cores):

```bash
nohup sh -c 'python -m benchmarks.liana_final --part all && python -m benchmarks.liana_final.amendment2 --part all' > run.log 2>&1 &
```

On a smaller machine, or split across people, run the parts separately. `--only REGEX` keeps only the runs whose tag
matches, so a tissue or a replicate can go to another machine.

| Step | Command | Needs | RAM |
|---|---|---|---|
| 1 | `python -m benchmarks.liana_final --part run --dataset <name>` (once per data set: `tnbc_wu2021`, `allen_brain_atlas`, `tabula_sapiens_mammary`, `tabula_sapiens_liver`, `moffitt2018_hypothalamus`) | data | 8 GB |
| 1 | `python -m benchmarks.liana_final --part calibration --dataset tabula_sapiens_mammary` (and `tabula_sapiens_liver`) | data | 8 GB |
| 1 | `python -m benchmarks.liana_final --part power` | data | 8 GB |
| 1 | `python -m benchmarks.liana_final.amendment2 --part <p>` for `hard_power`, `recal`, `mouse`, `timing` | data | 8 GB |
| 1 | `python -m benchmarks.liana_final.amendment2 --part confirm` and `--part donor --only '^tabula_sapiens_lung'` (Lung) | data | 32 GB |
| 2 | `python -m benchmarks.liana_final --part biology` | all five `run` tables | 8 GB |
| 2 | `python -m benchmarks.liana_final.amendment2 --part donor --only '^(?!tabula_sapiens_lung)'` | `run` tables (TNBC, Mammary, Liver) | 8 GB |
| 2 | `python -m benchmarks.liana_final.amendment2 --part sensitivity` (Mammary at 5,000 permutations: 16 GB or more) | `run` tables (Mammary, hypothalamus) | 8–16 GB |
| 2 | `python -m benchmarks.liana_final.amendment2 --part cases` | all five `run` tables | 8 GB |
| 3 | `python -m benchmarks.liana_final.amendment2 --part biology` | step 2 `biology`, the donor tables | 8 GB |

"`run` tables" are the files `data/liana_final/<name>_wide.pkl` and `<name>_store.pkl`. Copy `data/liana_final/`
between machines when a later step runs elsewhere.

## 4. Collect

Send your `benchmarks/results/liana_final/` folder (JSON files) to whoever assembles the results. They run:

```bash
python -m benchmarks.liana_final.merge_results <your results folder> --machine "<who, which machine>"
python -m benchmarks.liana_final --part report
python -m benchmarks.liana_final.amendment2 --part report
python -m benchmarks.liana_final.figures
```

The merge stops if a file was produced under a different plan or amendment (the SHA-256 is in every file), or if
the same run came out differently on two machines.
