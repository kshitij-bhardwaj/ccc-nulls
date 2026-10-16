# Tool workspaces

One folder per tool. **Each folder has a single owner, who is the only person who edits it**
(assigned within the team). Shared code lives in `src/ccc_nulls/` and changes there go
through a reviewed pull request.

| Folder | Tool | LIANA method | Statistic (Open Problems ranks by) | Native null (Null A)? | Native implementation |
|---|---|---|---|---|---|
| `cellphonedb/` | CellPhoneDB v5 | `cellphonedb` | `lr_means` | **Yes**: permutes cluster labels (1,000 by default) | Python package `cellphonedb` |
| `cellchat/` | CellChat | `cellchat` | `lr_probs` (native: `prob`) | **Yes**: permutes cell labels (100 by default) | CellChat 2.2.0.9001 ported to Python (`tools/cellchat/cellchat_py.py`, no R) |
| `natmi/` | NATMI | `natmi` | `spec_weight` (edge specificity) | No permutation p-value | Python scripts (`asrhou/NATMI`) |
| `connectome/` | Connectome | `connectome` | `scaled_weight` (`weight_sc`) | No permutation null on the edge score | R package `Connectome` (Seurat) |

The native-null column comes from each tool's paper and docs. Every owner should confirm it
against the version they install.

## The contract: two files per tool

`adapter.py` gets the tool's interactions into the common schema:

```python
TOOL = "cellphonedb"            # identifier used in results
LIANA_METHOD = "cellphonedb"    # liana-py method name
SCORE_NAME = "lr_means"         # statistic tested by the nulls

def liana_reference(adata, groupby, **kw) -> pd.DataFrame:   # working now: LIANA's version
def prepare_inputs(adata, groupby, out_dir) -> dict:         # dataset -> the tool's own input files
def run_native(inputs, out_dir, **kw) -> dict:               # run the original tool
def native_to_schema(native) -> pd.DataFrame:                # the tool's output -> common schema
```

`nulls.py` is this tool's own Null A and Null B, written by the owner on the tool's own score:

```python
def null_a(adata, interactions, groupby="label", *, B=1000, seed=0) -> pd.DataFrame:  # adds p_null_a
def null_b(adata, interactions, groupby="label", *, B=1000, seed=0) -> pd.DataFrame:  # adds p_null_b
```

- `interactions` is the common-schema table from `adapter.py`; `adata.X` is the expression the
  tool scored.
- Return the **same rows in the same order** with one p-value column added, every value in
  (0, 1] (untestable rows get 1). `ccc_nulls.schema.validate_null` checks this, in the tests
  and in the benchmark.
- Use a permutation p of the form `(1 + #{null >= observed}) / (B + 1)` so that p > 0.
- Shared helpers you may call (optional): `ccc_nulls.multitool.StatsEngine` and
  `permute_codes`, `ccc_nulls.nulls.gene_strata`, `ccc_nulls.stats.bh`. If you need something
  every tool could use, propose it for `src/ccc_nulls/` in a PR rather than copying it.

The pipeline is:

```
dataset (h5ad) --prepare_inputs--> native input files --run_native--> native output
      |                                                                     |
      +--liana_reference--> common schema <------- native_to_schema --------+
                                  |
                  nulls.null_a / nulls.null_b   (this folder, this tool's score)
                                  |
          benchmarks/run_tool.py: checks, calls, native vs ours, Open Problems metrics (shared)
```

**Common schema** (`ccc_nulls.schema`): `ligand, receptor, source, target, score, score_name,
p_native, ligand_complex, receptor_complex, tool`. One row per interaction, with the cell-type
labels exactly as in the dataset. Run `ccc_nulls.schema.validate(df)` before returning.

## Running and comparing

```bash
pytest tools/<tool>                                              # contract tests (skip until implemented)
python benchmarks/run_tool.py --tool <tool> --quick --B 50       # 8,000-cell smoke test
python benchmarks/run_tool.py --tool <tool> --B 1000             # full TNBC, LIANA interactions
python benchmarks/run_tool.py --tool <tool> --source native      # once the native adapter works
```

`run_tool.py` writes `benchmarks/results/tools/<tool>_<source>_<tag>.json` with, for every
p-value present (Null A, Null B, native): calls at p < 0.05 and BH q < 0.05; agreement of
ours with the native p (Spearman, calls in both / only one); and the Open Problems AUPRC and
odds ratio when the score is gated at p <= 0.05. Every tool is scored by the same code.

## Definition of done for a tool

1. `null_a` and `null_b` pass `pytest tools/<tool>` and run on TNBC through `run_tool.py`.
2. CellPhoneDB and CellChat: `p_null_a` agrees with the tool's own permutation p-value
   (`vs_native` in the summary). Explain any gap. This is plan step 3.
3. `native_to_schema(run_native(prepare_inputs(adata, ...)))` passes `schema.validate`, and is
   compared with `liana_reference` (overlap of interactions, score correlation). LIANA
   re-implements these tools, so some difference is expected; document it.

## Datasets

See the root `README.md` and `data/README.md`. Load the evaluation data with
`ccc_nulls.datasets.load_openproblems("tnbc_wu2021")`.
