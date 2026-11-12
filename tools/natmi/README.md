# NATMI

Owner: see `.github/CODEOWNERS`. Only the owner edits this folder; anything shared goes in
`src/ccc_nulls/` through a reviewed pull request.

| | |
|---|---|
| LIANA method | `natmi` (resource `connectomedb2020`) |
| Statistic tested by the nulls | `spec_weight`: NATMI's *edge average expression derived specificity* |
| Native null (Null A) | None: NATMI reports scores, not p-values. Null A here is the label permutation applied to its score |
| Native implementation | `ExtractEdges.py` of `github.com/asrhou/NATMI`, commit `f35f677` (2021-04-13), database `lrc2p` |
| Benchmark | `benchmarks/liana_nulls/AMENDMENT_5.md`, run by `python -m benchmarks.liana_nulls.amendment5` |

## The score

For ligand L in sender s and receptor R in receiver t, with means over all cells of a type (zeros included):

    spec_weight = (mean_L[s] / sum_k mean_L[k]) * (mean_R[t] / sum_k mean_R[k])
    expr_prod   = mean_L[s] * mean_R[t]

The ligand factor reads only the ligand's group means and the receptor factor only the receptor's, so the score
splits exactly into a sender and a receiver side (`adapter.side_score`), and Null B runs on it through
`ccc_nulls.composite`. lrc2p has no heteromeric complexes.

## Status

| Function in `adapter.py` | State |
|---|---|
| `liana_reference` | working (LIANA's re-implementation, on `connectomedb2020`) |
| `prepare_inputs` | working: `em.csv` (genes x cells, the database's genes), `ann.csv` |
| `run_native` | working; needs NATMI's own environment (below) |
| `native_to_schema` | working: `Edges_lrc2p.csv` -> common schema |
| `run`, `side_score`, `null_a`, `null_b`, `testable` | working: every (pair, source, target), Null A and Null B on one permutation sequence |
| `fidelity` | column 1 vs column 2: same positive rows, scores equal to 1e-9 (relative) |

Tested on toy data against NATMI's arithmetic transcribed from `ExtractEdges.py`. Not yet run against native
NATMI or on the benchmark's data sets.

## Native environment

NATMI needs the versions it was tested with (`requirements-native.txt`: pandas 1.0.3, xlrd 1.2.0); modern pandas
cannot run it.

```bash
python3.7 -m venv .venv-natmi
.venv-natmi/bin/pip install -r tools/natmi/requirements-native.txt
git clone https://github.com/asrhou/NATMI.git data/natmi/NATMI
git -C data/natmi/NATMI checkout f35f677cb2e4c8a176a6c5501a9483c9e28661a1
export NATMI_PYTHON=$PWD/.venv-natmi/bin/python NATMI_DIR=$PWD/data/natmi/NATMI
```

The Python side (`run`, the nulls) reads `lrdbs/lrc2p.csv` from the same clone.

## Run the tests

```bash
pytest tools/natmi
```

The native pipeline test skips when the clone is missing.
