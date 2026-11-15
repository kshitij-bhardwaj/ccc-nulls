PY ?= python
export PYTHONPATH := $(CURDIR)/src:$(CURDIR)

.PHONY: install install-r test data benchmark1 benchmark2 reports figures

install:            ## the locked Python environment that produced the published results
	$(PY) -m pip install -r requirements.txt && $(PY) -m pip install --no-deps -e .

install-r:          ## R + CellChat 2.1.2, only for Benchmark 1's native CellChat column
	Rscript tools/cellchat/install_r_packages.R

test:               ## no downloaded data needed; tests that need it skip themselves
	$(PY) -m pytest -q

data:               ## ~10 GB, md5-checked (TNBC and Allen brain need the copies described in data/README.md)
	$(PY) -m ccc_nulls.datasets tnbc_wu2021 allen_brain_atlas tabula_sapiens_mammary tabula_sapiens_liver tabula_sapiens_lung

benchmark1: data    ## Benchmark 1, each tool's own p-values (also needs R and the auxiliary files)
	$(PY) -m benchmarks.liana_nulls --part core
	$(PY) -m benchmarks.liana_nulls --part core --dataset moffitt2018_hypothalamus
	$(PY) -m benchmarks.liana_nulls --part calibration power biology amendment extend power_a3 b2_validity ablation report
	$(PY) -m benchmarks.liana_nulls.amendment4 --part all
	$(PY) -m benchmarks.liana_nulls.amendment5 --part all

benchmark2: data    ## Benchmark 2, LIANA+'s consensus (days on one 32 GB machine; see kaggle/README.md)
	$(PY) -m benchmarks.liana_final --part all
	$(PY) -m benchmarks.liana_final.amendment2 --part all
	$(PY) -m benchmarks.liana_final.amendment3 --part all
	$(PY) -m benchmarks.liana_final.amendment4 --part all
	$(PY) -m benchmarks.liana_final.amendment5 --part all
	$(PY) -m benchmarks.liana_final.amendment6 --part all
	$(PY) -m benchmarks.liana_final.amendment7 --part all
	$(PY) -m benchmarks.liana_final.amendment8 --part all

reports:            ## every results README (the tables) from the committed JSON files (no data, seconds)
	$(PY) -m benchmarks.liana_nulls --part report
	$(PY) -m benchmarks.liana_final --part report
	$(PY) -m benchmarks.liana_final.amendment2 --part report
	$(PY) -m benchmarks.liana_final.amendment3 --part report
	$(PY) -m benchmarks.liana_final.amendment4 --part report
	$(PY) -m benchmarks.liana_final.amendment6 --part report
	$(PY) -m benchmarks.liana_final.report_amendment7
	$(PY) -m benchmarks.liana_final.amendment8 --part report

figures:            ## the report's three figures from the committed JSON files (no data, seconds)
	$(PY) figures/make_figures.py
