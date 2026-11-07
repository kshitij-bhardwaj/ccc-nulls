# Tool workspaces

One folder per tool. **Each folder has a single owner, who is the only person who edits it**
(see `.github/CODEOWNERS`). Shared code lives in `src/ccc_nulls/` and changes there go
through a reviewed pull request.

| Folder | Tool | LIANA method | Statistic (Open Problems ranks by) | Native null (Null A)? | Native implementation |
|---|---|---|---|---|---|
| `cellphonedb/` | CellPhoneDB v5 | `cellphonedb` | `lr_means` | **Yes**: permutes cluster labels (1,000 by default) | Python package `cellphonedb` |
| `cellchat/` | CellChat | `cellchat` | `lr_probs` (native: `prob`) | **Yes**: permutes cell labels (100 by default) | CellChat 2.1.2 ported to Python (`cellchat/cellchat_py.py`), no R |
| `natmi/` | NATMI | `natmi` | `spec_weight` (edge specificity) | No permutation p-value | Python scripts (`asrhou/NATMI`) |
| `connectome/` | Connectome | `connectome` | `scaled_weight` (`weight_sc`) | No permutation null on the edge score | R package `Connectome` (Seurat) |

The native-null column comes from each tool's paper and docs. Every owner should confirm it
against the version they install.

## The contract: every `adapter.py` exposes the same five things

```python
TOOL = "cellphonedb"            # identifier used in results
LIANA_METHOD = "cellphonedb"    # liana-py method name
SCORE_NAME = "lr_means"         # statistic tested by the nulls

def liana_reference(adata, groupby, **kw) -> pd.DataFrame:   # working now: LIANA's version
def prepare_inputs(adata, groupby, out_dir) -> dict:         # dataset -> the tool's own input files
def run_native(inputs, out_dir, **kw) -> dict:               # run the original tool
def native_to_schema(native) -> pd.DataFrame:                # the tool's output -> common schema
```

The pipeline is:

```
dataset (h5ad) --prepare_inputs--> native input files --run_native--> native output
      |                                                                     |
      +--liana_reference--> common schema <------- native_to_schema --------+
                                  |
                      null models + evaluation (shared, src/ccc_nulls)
```

**Common schema** (`ccc_nulls.schema`): `ligand, receptor, source, target, score, score_name,
p_native, ligand_complex, receptor_complex, tool`. One row per interaction, with the cell-type
labels exactly as in the dataset. Run `ccc_nulls.schema.validate(df)` before returning.

### Optional: Null B for your tool's exact score

Null B (`ccc_nulls.composite`, the composite null "sender side specific AND receiver side
specific") needs one more function in the adapter: the tool's score split into its two sides.

```python
from ccc_nulls import composite

def side_score(adata, interactions, groupby, **kw) -> dict:
    side = composite.SideScore(
        summarize=lambda codes, K: ...,            # per-group statistics of the genes, for labels `codes`
        score=lambda sender, receiver: ...,        # the tool's score: sender side from `sender`,
    )                                              #   receiver side from `receiver`
    return dict(side=side, codes=..., levels=..., names=...)
```

`score(summarize(c), summarize(c))` must equal the tool's score. Then
`composite.get(method)(side, codes, K, B=..., seed=...)` runs any registered Null B on it;
`tools/cellchat/adapter.py` (`side_score`, `null_b`) is the worked example. A score computed from
cluster means needs no code: `composite.mean_side_score` or `composite.null_b_table`. New Null B
methods go in `ccc_nulls.composite` with `@composite.register("<name>")` and then run on every tool.

## Definition of done for a tool

1. `native_to_schema(run_native(prepare_inputs(adata, ...)))` passes `schema.validate`.
2. Compared with `liana_reference` on the same data: report the overlap of interactions and the
   score correlation, and explain any difference. LIANA re-implements these tools, so some
   difference is expected.
3. Tools with a native null (CellPhoneDB, CellChat): the native p-values are compared with
   our Null A (`ccc_nulls.multitool.label_null`). This is plan step 3.
4. `pytest tools/<tool>` passes.

## Running the pipeline

The settings below are frozen: they were chosen on Tabula Sapiens and then applied once to TNBC.
CellChat is the worked example; the exact numbers are in `cellchat/README.md`, not repeated here.

### Environment

- Run from the repository root with the repository's environment, as modules:
  `.venv/bin/python -m tools.cellchat.<script>`, `.venv/bin/python -m pytest`. The scripts import
  `tools.<tool>` and `ccc_nulls`, so running a file by its path from another folder fails.
- No R. CellChat is the Python port in `cellchat/cellchat_py.py`; `pip install rdata` once, to read
  CellChatDB from the R package's data file.
- Data go to `data/` (or `$CCC_NULLS_DATA`):

```bash
.venv/bin/python -m ccc_nulls.datasets tnbc_wu2021 tabula_sapiens_mammary tabula_sapiens_liver
```

### Inputs

- **Expression:** `adata.X` is log-normalised, ln(1 + CP10k), cells × genes, not scaled.
  `datasets.load_openproblems("tnbc_wu2021")` loads TNBC's `normalized` layer;
  `datasets.load_tabula_sapiens("tabula_sapiens_<tissue>")` normalises Tabula Sapiens' raw counts the
  same way. Both drop cell types with fewer than 5 cells.
- **Groups:** `groupby="label"` for both (Tabula Sapiens' `label` is its `cell_type`).
- **CellChat version:** 2.1.2, the release, is the default. Null B's side score exists for 2.x only.
- **Candidate rows:** `adapter.native_to_schema(adapter.run_native(adapter.prepare_inputs(adata, groupby), nboot=2, seed=1))`.
  This gives every L-R pair × sender × receiver with CellChat's `prob`. `nboot=2` only builds the
  rows: its `p_native` comes from 2 permutations and is not a p-value to use. The nulls compute p.

### Null settings

| Setting | Value |
|---|---|
| Permutations | B = 1,000 |
| Same permutations for both nulls | the same `B`, `seed` and `rng` in `nulls.null_a(...)` and `adapter.null_b(...)`. Pass `rng` explicitly: `null_a` defaults to `"cellchat"`, `null_b` to `"numpy"` |
| Null B | `method="composite"`, `delta=0` (δ > 0 is a sensitivity run only) |
| Ties | the tool's own convention, for both nulls. CellChat: p = (1 + #{s_b > s}) / (B + 1), strict `>`; `adapter.null_b` passes `ties="ignore"`. `composite`'s default `ties="count"` (≥) is for tools whose own test counts ties |
| Multiple testing | BH over `adapter.testable(null_a_result, null_b_result)` rows only: some labelling gives score > 0. Not over every row, and not over the rows the tool reports (score > 0 with the observed labels), which select on the observed labels |
| Level | α = 0.05. α = 0.01 and Benjamini–Yekutieli need B ≥ 5,000: the smallest p is 1/(B + 1) |

### Calibration (Tabula Sapiens, known-null data)

```bash
.venv/bin/python -m tools.cellchat.calibrate_null_a --tissue all                                   # Null A, labels fully shuffled
.venv/bin/python -m tools.cellchat.null_b_calibration --tissue mammary --scenario receptor --reps 20
.venv/bin/python -m tools.cellchat.null_b_calibration --tissue mammary --scenario ligand --reps 20  # also labels, real; --tissue liver
.venv/bin/python -m tools.cellchat.null_b_calibration --summarise
.venv/bin/python -m tools.cellchat.null_b_checklist --part pcal
```

- **Scenarios:**
  - `labels`: every cell's label shuffled, so both sides are null.
  - `receptor`: the receptor-side genes (receptors, co-receptors) are shuffled across cells with one
    permutation; the ligand side keeps its cell-type signal.
  - `ligand`: the mirror image.
  - `real`: the observed labels, for the number of calls.
- **Known-null rows:** in `receptor`, every row of an L-R pair whose receiver-side genes were all
  shuffled, in every (sender, receiver); `ligand` likewise. Genes on both sides (agonists,
  antagonists, a gene that is a ligand in one pair and a receptor in another) are never shuffled, so
  pairs that use them are not known null; they stay in the data and in the BH family.
- **FDR:** the mean over replicates of V / R, where V is calls on known-null rows and R is all calls.
  It is a lower bound, since calls on the other rows can be false too. On `labels` every call is
  false, so the measure there is P(any call).

### TNBC is held out

Freeze every setting on Tabula Sapiens first, then run TNBC once:
`.venv/bin/python -m tools.cellchat.null_b_tnbc`. Do not change a setting after seeing TNBC results.
A change means calibrating again on Tabula Sapiens and a new, documented TNBC run.

### Pitfalls

- **Two functions are called Null B.** Null B is `ccc_nulls.composite` (`method="composite"`, through
  the adapter's `null_b`). `ccc_nulls.nulls.null_b` is the earlier abundance-matched decoy null, kept
  for the record, and is not Null B.
- **B2 is that decoy null with its premise repaired** (`ccc_nulls.competitive`, Amendment 3 of
  `benchmarks/liana_nulls`): the inter-gene correlation goes into the variance, as in CAMERA. It
  needs no side score, works on any tool's rows, and takes the decoys from genes outside every
  database, so data loaded as database genes only (Tabula Sapiens) needs extra genes loaded.
- **`SideScore` requirements:**
  - `score(summarize(c, K), summarize(c, K))` must equal the tool's score exactly;
  - the score is non-negative, with 0 meaning no signal (p = 1 there);
  - `summarize(codes, K)` takes 0-based group codes;
  - the sender side is read only from `score`'s first argument and the receiver side only from its second.
- **Permutations are 0-based in `composite`** (a cells × B matrix of cell indices). CellChat's
  R-style permutations (`cellchat_py.r_permutations`, `nulls.null_a(permutations=...)`) are 1-based;
  `adapter.null_b` converts them. Do not pass one to the other unconverted.
- **Do not use B1** (`nulls.null_b_empirical`, Efron's empirical null). It is not identifiable here,
  because most rows depend on cell type; it is kept as a negative result (`cellchat/README.md`).

## Datasets

See the root `README.md` and `data/README.md`. Load the evaluation data with
`ccc_nulls.datasets.load_openproblems("tnbc_wu2021")`, and the calibration tissues with
`ccc_nulls.datasets.load_tabula_sapiens("tabula_sapiens_mammary")` (or `_liver`).
