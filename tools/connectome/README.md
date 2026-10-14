# Connectome

One owner, assigned within the team. Only the owner edits this folder; anything shared goes in
`src/ccc_nulls/` through a reviewed pull request.

| | |
|---|---|
| LIANA method | `connectome` |
| Statistic tested by the nulls | `scaled_weight` |
| Native null (Null A) | No permutation null on the edge score; it reports per-gene differential-expression p-values (ligand, receptor) |
| Native implementation | R package `Connectome` on a Seurat object (run from Python through `Rscript`) |

## Status

| Function in `adapter.py` | State |
|---|---|
| `liana_reference` | working (LIANA's re-implementation, common schema) |
| `prepare_inputs` | to do |
| `run_native` | to do |
| `native_to_schema` | to do |

| Function in `nulls.py` | State |
|---|---|
| `null_a` | to do: shuffle labels, recompute this tool's score |
| `null_b` | to do: abundance-matched decoy genes, this tool's score |

## Native input (to verify against the installed version)

1. Expression for a Seurat object: counts or log-normalised matrix + cell identities.
    2. Seurat preprocessing: `NormalizeData`, `ScaleData` (on the LR genes at least), `Idents()` = cell type.

## How the native tool is run (to verify)

`CreateConnectome(object, species = "human", p.values = TRUE, calculate.DOR = FALSE)`; write the edge table to CSV.

## Native output (to verify)

Edge table with `source, target, ligand, receptor, weight_norm, weight_sc, p_val_lig, p_val_rec, ...`.

## Checklist

- [ ] confirm column names and default filters in the installed version
- [ ] use `weight_sc` as `score` (= LIANA's scaled_weight)
- [ ] decide how to report the DE p-values (they are not a null for the edge score)
- [ ] `nulls.null_a` (our implementation; no native p to compare)
- [ ] `nulls.null_b`

## Run the tests and the benchmark

```bash
pytest tools/connectome
python benchmarks/run_tool.py --tool connectome --quick --B 50     # smoke test, 8,000 cells
python benchmarks/run_tool.py --tool connectome --B 1000             # full TNBC run
```
