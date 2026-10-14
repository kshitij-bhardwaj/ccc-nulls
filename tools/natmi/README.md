# NATMI

One owner, assigned within the team. Only the owner edits this folder; anything shared goes in
`src/ccc_nulls/` through a reviewed pull request.

| | |
|---|---|
| LIANA method | `natmi` |
| Statistic tested by the nulls | `spec_weight` |
| Native null (Null A) | None: natmi reports expression weights and specificities, not p-values |
| Native implementation | Python scripts from `github.com/asrhou/NATMI` (`ExtractEdges.py`) |

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

1. `em.csv`: log-normalised expression, genes x cells, first column = gene symbol.
    2. `ann.csv`: two columns, cell barcode and cell-type annotation.
    3. Interaction database: `lrc2p` (connectomeDB2020; ships with NATMI).

## How the native tool is run (to verify)

`python ExtractEdges.py --species human --emFile em.csv --annFile ann.csv --interDB lrc2p --out <dir>`

## Native output (to verify)

`Edges_lrc2p.csv`: `Sending cluster, Ligand symbol, Receptor symbol, Target cluster`,
    detection rates, `Edge average expression weight`, `Edge average expression derived specificity`, ...

## Checklist

- [ ] confirm column names in the installed NATMI version
- [ ] use `Edge average expression derived specificity` as `score` (= LIANA's spec_weight)
- [ ] no native null: Null A here is our implementation only; say so in results
- [ ] `nulls.null_a` (our implementation; no native p to compare)
- [ ] `nulls.null_b`

## Run the tests and the benchmark

```bash
pytest tools/natmi
python benchmarks/run_tool.py --tool natmi --quick --B 50     # smoke test, 8,000 cells
python benchmarks/run_tool.py --tool natmi --B 1000             # full TNBC run
```
