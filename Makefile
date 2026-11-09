PY ?= python

.PHONY: install test docs data reproduce

install:            ## the environment that produced the published results
	$(PY) -m pip install -r requirements/lock-liana-final.txt && $(PY) -m pip install -e . --no-deps

test:
	$(PY) -m pytest -q --ignore=tools/connectome

docs:
	sphinx-build -b html docs docs/_build/html

data:               ## ~9 GB, md5-checked (the hypothalamus: see data/README.md)
	$(PY) -m ccc_nulls.datasets tnbc_wu2021 allen_brain_atlas tabula_sapiens_mammary tabula_sapiens_liver tabula_sapiens_lung

reproduce: data     ## every result in benchmarks/results/liana_final (about a day and a half on 8 GB)
	$(PY) -m benchmarks.liana_final --part all
	$(PY) -m benchmarks.liana_final.amendment2 --part all
