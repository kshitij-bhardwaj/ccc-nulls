# Tool workspaces

One folder per tool. **Each folder has a single owner, who is the only person who edits it**
(see `.github/CODEOWNERS`). Shared code lives in `src/ccc_nulls/` and changes there go
through a reviewed pull request.

| Folder | Tool | LIANA method | Statistic (Open Problems ranks by) | Native null (Null A)? | Native implementation |
|---|---|---|---|---|---|
| `cellphonedb/` | CellPhoneDB v5 | `cellphonedb` | `lr_means` | **Yes**: permutes cluster labels (1,000 by default) | Python package `cellphonedb` |
| `cellchat/` | CellChat | `cellchat` | `lr_probs` | **Yes**: permutes group labels (100 by default) | R package `CellChat` |
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

## Definition of done for a tool

1. `native_to_schema(run_native(prepare_inputs(adata, ...)))` passes `schema.validate`.
2. Compared with `liana_reference` on the same data: report the overlap of interactions and the
   score correlation, and explain any difference. LIANA re-implements these tools, so some
   difference is expected.
3. Tools with a native null (CellPhoneDB, CellChat): the native p-values are compared with
   our Null A (`ccc_nulls.multitool.label_null`). This is plan step 3.
4. `pytest tools/<tool>` passes.

## Datasets

See the root `README.md` and `data/README.md`. Load the evaluation data with
`ccc_nulls.datasets.load_openproblems("tnbc_wu2021")`.
