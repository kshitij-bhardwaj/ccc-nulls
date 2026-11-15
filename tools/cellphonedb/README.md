# CellPhoneDB (v5)

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
| `prepare_inputs` | sparse float32 `counts.h5ad` + `meta.tsv`; optional subset to database genes |
| `run_native` | Method 1 (`cpdb_analysis_method.call`), means only |
| `native_to_schema` | Method 1 means -> common schema (`p_native` is NaN) |
| `method1_with_nulls` / `combine` | Method 1 + this project's Null A on CellPhoneDB's own statistic |

`NATIVE_MODE = "method1"`: the native path produces no native p-values. `HAS_NATIVE_NULL`
stays `True` because CellPhoneDB itself has one (Method 2), which is not run here.

## Native input

1. `counts.h5ad`: normalised, unscaled expression (cells x genes, HGNC symbols as `var_names`),
   sparse float32. CellPhoneDB reads `.h5ad` directly.
2. `meta.tsv`: two columns `Cell`, `cell_type`.
3. The CellPhoneDB database zip, downloaded with
   `cellphonedb.utils.db_utils.download_database(out_dir, "v5.0.0")`.

Cell IDs may contain dashes. Cell-type labels must not contain `|`, CellPhoneDB's
pair-column separator. Pass `database_path=` to `prepare_inputs` to write only genes the
database knows. CellPhoneDB drops all other genes anyway, but it loads the counts densely,
so on a full atlas this is the difference between ~5 GB and a few hundred MB.

## The statistic (matches cellphonedb 5.0.1)

`cpdb_null.py` reproduces CellPhoneDB's Method 1 / Method 2 statistic. Its `load()` reuses
CellPhoneDB's own loading and prefiltering, so it sees exactly CellPhoneDB's data.

- Counts are mapped to CellPhoneDB multidata ids as float32. Identifiers that map to the
  same id are averaged and all-zero genes are dropped. An interaction is kept only if every
  partner, including every subunit of a complex, is present.
- For each cell type, the mean and the fraction of expressing cells (> 0) are computed.
  **A complex takes the minimum over its subunits**, for both.
- Score of partner A in type x and partner B in type y: `(m_A > 0) * (m_B > 0) * (m_A + m_B) / 2`.
- Expression gate: `frac_A > threshold and frac_B > threshold` (strict, default 0.1).

`test_cpdb_null.py::test_matches_cellphonedb_helpers` checks means, fractions, scores and
the gate against CellPhoneDB's own `build_clusters`, `mean_analysis` and `percent_analysis`.
`combine()` also checks every recomputed score against CellPhoneDB's `means` output
(within its 3-decimal rounding) and raises if any differ.

## Null A on that statistic

```python
from tools.cellphonedb.adapter import method1_with_nulls

res = method1_with_nulls(
    adata, "cell_type", "outputs/cellphonedb",
    database_path="/path/to/cellphonedb.zip", B=1000, threshold=0.1, seed=0,
)
res.loc[res.tested, ["ligand", "receptor", "source", "target", "score", "p_rand", "q_rand"]]
```

- Cell-type labels are permuted uniformly (our own NumPy generator), and the same score
  (complex minimum included) is recomputed for every permutation. As in CellPhoneDB
  Method 2, permuted scores are not gated.
- `p_rand = (1 + #{b : s_b >= s_obs}) / (B + 1)`. `pvalue="cpdb"` gives CellPhoneDB's
  `#{b : s_b > s_obs} / B` instead, for side-by-side checks.
- Only **tested** rows (gate passed, score > 0) get a p-value. Benjamini–Hochberg runs over
  tested rows only. Other rows have `p_rand = q_rand = NaN`. CellPhoneDB would give them p = 1,
  but they are left out here so they do not inflate the number of tests.
- One row per CellPhoneDB interaction x (source, target). For a complex, `ligand` /
  `receptor` is the subunit that carries the complex's minimum mean in that cell type, and
  `ligand_complex` / `receptor_complex` is the complex name.
- Interactions where both or neither partner is a receptor have no ligand -> receptor
  direction and are skipped, with a warning that gives the count.
- Null B is not computed on this path yet.

## Native output (verified against 5.0.1 source)

`means_result`: one row per interaction (`id_cp_interaction`, `interacting_pair`,
`partner_a`, `partner_b`, `gene_a`, `gene_b`, `receptor_a`, `receptor_b`, ...) and one column
per `"<type a>|<type b>"`, rounded to 3 decimals. Partner A is in type a, partner B in type b.
Partner order is **not** ligand-first, so `receptor_a` / `receptor_b` decide the direction.
`gene_a` / `gene_b` are empty for complexes. The `deconvoluted` table has **no** `receptor`
column, and its means are complex-level, so it cannot be used to pick subunits. `cpdb_null`
computes them instead.

## Checklist

- [x] confirm v5 input/output format against the installed version (5.0.1 source)
- [x] orient partners (ligand vs receptor) using receptor_a / receptor_b
- [x] complexes: minimum-subunit rule; record the subunit used
- [x] statistic identical to CellPhoneDB's (tested against its helpers)
- [ ] run `method1_with_nulls` end to end with a downloaded database (`CPDB_DB=... pytest tools/cellphonedb`)
- [ ] optional: compare `pvalue="cpdb"` Null A with Method 2 `pvalues`

## Run the tests

```bash
pytest tools/cellphonedb                                 # no database needed
CPDB_DB=/path/to/cellphonedb.zip pytest tools/cellphonedb  # + native CellPhoneDB runs
```
