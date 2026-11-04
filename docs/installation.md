# Installation

ccc-nulls needs Python 3.10 or newer. The core library depends on numpy, scipy and pandas only; the
LIANA+ integration (`ccc_nulls.rank_aggregate`) needs the `liana` extra.

```bash
pip install "ccc-nulls[liana] @ git+https://github.com/kshitij-bhardwaj/ccc-nulls.git"
```

LIANA+ 1.8 needs anndata < 0.12, zarr < 3, pandas < 3 and numpy < 2.3; the extra pins these.

## Exact environment of the published results

`requirements/lock-liana-final.txt` is the environment every number in the benchmark was produced in
(Python 3.12, LIANA+ 1.8.0, anndata 0.11.4, scanpy 1.11.5):

```bash
git clone https://github.com/kshitij-bhardwaj/ccc-nulls.git && cd ccc-nulls
python -m venv .venv && source .venv/bin/activate
pip install -r requirements/lock-liana-final.txt && pip install -e . --no-deps
pytest -q --ignore=tools/connectome
```

## Development

```bash
pip install -e ".[liana,test,docs]"
pytest -q --ignore=tools/connectome
sphinx-build -b html docs docs/_build/html
```
