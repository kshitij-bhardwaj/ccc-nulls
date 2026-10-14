# CellPhoneDB (v5)

Owner: see `.github/CODEOWNERS`. Only the owner edits this folder; anything shared goes in
`src/ccc_nulls/` through a reviewed pull request.

| | |
|---|---|
| LIANA method | `cellphonedb` |
| Statistic tested by the nulls | `lr_means` |
| Native null (Null A) | Yes, permutes cluster labels (default 1,000 iterations); reported as `pvalues` |
| Native implementation | Python package `cellphonedb` (pip install cellphonedb) |

## Status

| Function in `adapter.py` | State |
|---|---|
| `liana_reference` | working (LIANA's re-implementation, common schema) |
| `prepare_inputs` | to do |
| `run_native` | to do |
| `native_to_schema` | to do |

## Native input (to verify against the installed version)

1. `counts.h5ad`: log-normalised expression, cells x genes, gene symbols as var_names.
    2. `meta.tsv`: two columns `barcode_sample`, `cell_type` (cell barcode -> label).
    3. The CellPhoneDB database zip (`cellphonedb.zip`), downloaded with
       `cellphonedb.utils.db_utils.download_database(out_dir, "v5.0.0")`.

## How the native tool is run (to verify)

`cpdb_statistical_analysis_method.call(cpdb_file_path=..., meta_file_path=..., counts_file_path=...,
    counts_data="hgnc_symbol", output_path=..., iterations=1000, threshold=0.1)`

## Native output (to verify)

Wide tables `means`, `pvalues`, `significant_means`: one row per interacting pair
    (`interacting_pair`, `partner_a`, `partner_b`, `gene_a`, `gene_b`, `receptor_a`, `receptor_b`, ...)
    and one column per `"<source>|<target>"`. Partner order is **not** ligand-first: use the
    `receptor_a` / `receptor_b` flags to decide which partner is the ligand.

## Checklist

- [ ] confirm v5 input/output format against the installed version
- [ ] orient partners (ligand vs receptor) using receptor_a / receptor_b
- [ ] complexes: map `complex:...` partners to their subunits; record the subunit used
- [ ] compare native p-values with Null A (`ccc_nulls.multitool.label_null`)

## Run the tests

```bash
pytest tools/cellphonedb
```
