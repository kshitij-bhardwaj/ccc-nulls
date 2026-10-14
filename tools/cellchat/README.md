# CellChat

Owner: see `.github/CODEOWNERS`. Only the owner edits this folder; anything shared goes in
`src/ccc_nulls/` through a reviewed pull request.

| | |
|---|---|
| LIANA method | `cellchat` |
| Statistic tested by the nulls | `lr_probs` |
| Native null (Null A) | Yes, permutes group labels (`nboot = 100` by default in `computeCommunProb`); reported as `pval` |
| Native implementation | R package `CellChat` (run from Python through `Rscript`) |

## Status

| Function in `adapter.py` | State |
|---|---|
| `liana_reference` | working (LIANA's re-implementation, common schema) |
| `prepare_inputs` | to do |
| `run_native` | to do |
| `native_to_schema` | to do |

## Native input (to verify against the installed version)

1. Expression: log-normalised genes x cells, written as `matrix.mtx` + `genes.tsv` + `barcodes.tsv`.
    2. `meta.csv`: cell barcode -> `labels`.
    3. Database: `CellChatDB.human` (ships with the R package).

## How the native tool is run (to verify)

An R script (e.g. `run_cellchat.R`): `createCellChat(object = data, meta = meta, group.by = "labels")`,
    `subsetData`, `identifyOverExpressedGenes`, `identifyOverExpressedInteractions`,
    `computeCommunProb(type = "triMean")`, `filterCommunication(min.cells = 10)`, then
    `subsetCommunication(thresh = 1)` to keep non-significant rows too; write it to CSV.

## Native output (to verify)

Long table: `source, target, ligand, receptor, prob, pval, interaction_name, ...`.

## Checklist

- [ ] pin the CellChat version and record it in this README
- [ ] keep all rows (thresh = 1), not only significant ones, so nulls can be compared
- [ ] check that `ligand` / `receptor` names match gene symbols (complexes use underscores)
- [ ] compare native `pval` with Null A (`ccc_nulls.multitool.label_null`)

## Run the tests

```bash
pytest tools/cellchat
```
