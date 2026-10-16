# Contributing

## Ownership

| Path | Who edits | How changes land |
|---|---|---|
| `tools/<tool>/` | that tool's owner only (assigned within the team) | PR, owner merges |
| `src/ccc_nulls/`, `benchmarks/`, `tests/` | anyone | PR + review by one other member |
| `README.md`, `docs/`, `CONTRIBUTING.md` | anyone | PR + review by one other member |

If your tool needs something in the shared library (a new helper, a schema field), open a PR
against `src/ccc_nulls/` instead of copying code into your folder.

## Branches

- `main` is always green: `pytest` passes.
- Work on a branch named `<area>/<short-description>`:
  `cellphonedb/native-adapter`, `cellchat/null-b`, `core/standardised-decoys`,
  `docs/tnbc-results`.
- Rebase on `main` before opening a PR; keep PRs small and about one thing.

## Commit messages

[Conventional Commits](https://www.conventionalcommits.org/), with the area as scope:

```
feat(cellphonedb): write counts and meta files in v5 format
fix(core): exclude candidate genes from the decoy pool
test(natmi): cover edge-table parsing
docs(readme): add spatial dataset link
refactor(evaluation): share Open Problems metrics across benchmarks
chore(deps): pin anndata below 0.12 for liana 1.8
```

Types: `feat`, `fix`, `test`, `docs`, `refactor`, `perf`, `chore`. Subject in the imperative,
lower case, no trailing period, at most ~72 characters. Use the body to explain *why*.

## Before you push

```bash
pytest                      # everything
pytest tools/<your-tool>    # just your workspace
```

- Never commit data, large results (`rows_*.csv.gz`, `*.pkl`) or notebooks with outputs.
  `data/` is git-ignored.
- Results that support a slide or report go in `benchmarks/results/` as small JSON or CSV,
  together with the command that produced them.
- No personal names in code, file names or commit messages.

## Environment

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[liana,cellchat,test]"
python -m ccc_nulls.datasets tnbc_wu2021
```

`liana` 1.8 needs `anndata<0.12`, `zarr<3`, `pandas<3`, `numpy<2.3`, `scanpy<1.12`; the
`[liana]` extra pins these. R-based tools (CellChat, Connectome) also need R ≥ 4.3 with
the package installed. Note the R and package versions in your tool's README.
