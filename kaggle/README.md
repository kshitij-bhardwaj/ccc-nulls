# The benchmark on Kaggle

Every result in `benchmarks/results/liana_final/` was produced by one of the Kaggle sessions below (CPU, 4 cores, 32 GB), from one commit of this repository, in the locked environment (`requirements/lock-liana-final.txt`). Each result file, and each run inside AMENDMENT_2's files, records the link of the session that produced it (`machine`). To rerun: `python kaggle/kernels.py code && python kaggle/kernels.py build && python kaggle/kernels.py dispatch && python kaggle/kernels.py collect`. To run elsewhere: [`RUN.md`](../RUN.md).

Inputs:

- code: [ccc-nulls-code](https://www.kaggle.com/datasets/kshitijbhardwaj1729/ccc-nulls-code) (`git archive` of the commit in `COMMIT`)
- auxiliary files: [ccc-nulls-aux](https://www.kaggle.com/datasets/kshitijbhardwaj1729/ccc-nulls-aux) (md5: `benchmarks/liana_final/aux_files.md5`)
- h5ad files: the output of [ccc-nulls-data](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-data), downloaded from Open Problems and CELLxGENE, md5-checked

| Session | Needs | Command | Run time |
|---|---|---|---|
| [data](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-data) | – | download the five h5ad files |  |
| [run-tnbc](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-run-tnbc) | – | `python -m benchmarks.liana_final --part run --dataset tnbc_wu2021` |  |
| [run-allen](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-run-allen) | – | `python -m benchmarks.liana_final --part run --dataset allen_brain_atlas` |  |
| [run-mammary](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-run-mammary) | – | `python -m benchmarks.liana_final --part run --dataset tabula_sapiens_mammary` |  |
| [run-liver](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-run-liver) | – | `python -m benchmarks.liana_final --part run --dataset tabula_sapiens_liver` |  |
| [run-hypothalamus](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-run-hypothalamus) | – | `python -m benchmarks.liana_final --part run --dataset moffitt2018_hypothalamus` |  |
| [calibration-mammary](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-calibration-mammary) | – | `python -m benchmarks.liana_final --part calibration --dataset tabula_sapiens_mammary` |  |
| [calibration-liver](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-calibration-liver) | – | `python -m benchmarks.liana_final --part calibration --dataset tabula_sapiens_liver` |  |
| [power](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-power) | – | `python -m benchmarks.liana_final --part power` |  |
| [biology](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-biology) | run-tnbc, run-allen, run-mammary, run-liver, run-hypothalamus | `python -m benchmarks.liana_final --part biology` |  |
| [a2-lung-receptor-0](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-receptor-0) | – | `python -m benchmarks.liana_final.amendment2 --part confirm --only '^tabula_sapiens_lung\\|receptor\\|0$'` |  |
| [a2-lung-receptor-1](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-receptor-1) | – | `python -m benchmarks.liana_final.amendment2 --part confirm --only '^tabula_sapiens_lung\\|receptor\\|1$'` |  |
| [a2-lung-receptor-2](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-receptor-2) | – | `python -m benchmarks.liana_final.amendment2 --part confirm --only '^tabula_sapiens_lung\\|receptor\\|2$'` |  |
| [a2-lung-ligand-0](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-ligand-0) | – | `python -m benchmarks.liana_final.amendment2 --part confirm --only '^tabula_sapiens_lung\\|ligand\\|0$'` |  |
| [a2-lung-ligand-1](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-ligand-1) | – | `python -m benchmarks.liana_final.amendment2 --part confirm --only '^tabula_sapiens_lung\\|ligand\\|1$'` |  |
| [a2-lung-ligand-2](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-ligand-2) | – | `python -m benchmarks.liana_final.amendment2 --part confirm --only '^tabula_sapiens_lung\\|ligand\\|2$'` |  |
| [a2-lung-power](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-power) | – | `python -m benchmarks.liana_final.amendment2 --part confirm --only '^tabula_sapiens_lung\\|f2\\|'` |  |
| [a2-lung-donor-global-0](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-global-0) | – | `python -m benchmarks.liana_final.amendment2 --part donor --only '^tabula_sapiens_lung\\|global\\|0$'` |  |
| [a2-lung-donor-global-1](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-global-1) | – | `python -m benchmarks.liana_final.amendment2 --part donor --only '^tabula_sapiens_lung\\|global\\|1$'` |  |
| [a2-lung-donor-global-2](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-global-2) | – | `python -m benchmarks.liana_final.amendment2 --part donor --only '^tabula_sapiens_lung\\|global\\|2$'` |  |
| [a2-lung-donor-donor-0](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-donor-0) | – | `python -m benchmarks.liana_final.amendment2 --part donor --only '^tabula_sapiens_lung\\|donor\\|0$'` |  |
| [a2-lung-donor-donor-1](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-donor-1) | – | `python -m benchmarks.liana_final.amendment2 --part donor --only '^tabula_sapiens_lung\\|donor\\|1$'` |  |
| [a2-lung-donor-donor-2](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-lung-donor-donor-2) | – | `python -m benchmarks.liana_final.amendment2 --part donor --only '^tabula_sapiens_lung\\|donor\\|2$'` |  |
| [a2-hard-power](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-hard-power) | – | `python -m benchmarks.liana_final.amendment2 --part hard_power` |  |
| [a2-recal-mammary](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-recal-mammary) | – | `python -m benchmarks.liana_final.amendment2 --part recal --only '^tabula_sapiens_mammary\\|'` |  |
| [a2-recal-liver](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-recal-liver) | – | `python -m benchmarks.liana_final.amendment2 --part recal --only '^tabula_sapiens_liver\\|'` |  |
| [a2-mouse](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-mouse) | – | `python -m benchmarks.liana_final.amendment2 --part mouse` |  |
| [a2-timing](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-timing) | – | `python -m benchmarks.liana_final.amendment2 --part timing` |  |
| [a2-donor-tnbc](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-donor-tnbc) | – | `python -m benchmarks.liana_final.amendment2 --part donor --only '^tnbc_wu2021\\|'` |  |
| [a2-donor-mammary](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-donor-mammary) | – | `python -m benchmarks.liana_final.amendment2 --part donor --only '^tabula_sapiens_mammary\\|'` |  |
| [a2-donor-liver](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-donor-liver) | – | `python -m benchmarks.liana_final.amendment2 --part donor --only '^tabula_sapiens_liver\\|'` |  |
| [a2-donor-real](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-donor-real) | run-tnbc, run-mammary, run-liver | `python -m benchmarks.liana_final.amendment2 --part donor --only '^(tnbc_wu2021\|tabula_sapiens_mammary\|tabula_sapiens_liver)$'` |  |
| [a2-sensitivity-seed-expr](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-sensitivity-seed-expr) | run-mammary | `python -m benchmarks.liana_final.amendment2 --part sensitivity --only '^seed\|@expr_prop='` |  |
| [a2-sensitivity-cells-labels](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-sensitivity-cells-labels) | run-mammary | `python -m benchmarks.liana_final.amendment2 --part sensitivity --only '@(cells\|labels)='` |  |
| [a2-sensitivity-perms-hypothalamus](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-sensitivity-perms-hypothalamus) | run-mammary, run-hypothalamus | `python -m benchmarks.liana_final.amendment2 --part sensitivity --only '^real$\|^moffitt2018_hypothalamus@n_perms'` |  |
| [a2-sensitivity-perms-mammary](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-sensitivity-perms-mammary) | run-mammary, run-hypothalamus | `python -m benchmarks.liana_final.amendment2 --part sensitivity --only '^real@\|^tabula_sapiens_mammary@n_perms'` |  |
| [a2-cases](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-cases) | run-tnbc, run-allen, run-mammary, run-liver, run-hypothalamus | `python -m benchmarks.liana_final.amendment2 --part cases` |  |
| [a2-biology](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-a2-biology) | biology, a2-donor-real, run-tnbc, run-allen, run-mammary, run-liver, run-hypothalamus | `python -m benchmarks.liana_final.amendment2 --part biology` |  |
| [report](https://www.kaggle.com/code/kshitijbhardwaj1729/ccc-nulls-report) | all | `python -m benchmarks.liana_final --part report`<br>`python -m benchmarks.liana_final.amendment2 --part report`<br>`python -m benchmarks.liana_final.figures` |  |
