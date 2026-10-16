# Data

Everything in this folder except this file is git-ignored. Download with:

```bash
python -m ccc_nulls.datasets tnbc_wu2021          # evaluation, cytokine truth (1.27 GB)
python -m ccc_nulls.datasets allen_brain_atlas    # evaluation, spatial truth (2.05 GB)
```

The downloader fetches 16 byte ranges in parallel and checks the md5. A single stream from S3
can be as slow as ~100 KB/s. A browser download of the link below works too; save the file
under the filename shown.

| Name | File | md5 | Link |
|---|---|---|---|
| `tnbc_wu2021` | `tnbc_wu2021_log_cp10k.h5ad` | `0c1c558c4137571eb5e9351c8408dc04` | https://openproblems-data.s3.amazonaws.com/resources/datasets/openproblems_v1/tnbc_wu2021/log_cp10k/dataset.h5ad |
| `allen_brain_atlas` | `allen_brain_atlas_log_cp10k.h5ad` | `89a5f0919ed3f2b3e03e7550ef5a20ef` | https://openproblems-data.s3.amazonaws.com/resources/datasets/openproblems_v1/allen_brain_atlas/log_cp10k/dataset.h5ad |
| `tabula_sapiens` | `tabula_sapiens_<tissue>.h5ad` | – | https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5 (one file per tissue) |
| Tabula Sapiens - Mammary (calibration tissue in use) | `tabula_sapiens_mammary.h5ad` | `6261c8dd5785c52fb8a18d0974602597` | https://datasets.cellxgene.cziscience.com/33b25f68-a9d6-4643-9283-653e937ae608.h5ad (1.33 GB, 30,936 cells; raw counts in `raw/X`, symbols in `raw/var/feature_name`, labels in `obs/cell_type`) |

## What is inside the Open Problems files

| Slot | TNBC (`tnbc_wu2021`) | Brain (`allen_brain_atlas`) |
|---|---|---|
| cells × genes | 42,512 × 28,200 (human) | 14,249 × 34,617 (mouse) |
| `layers["counts"]` | raw counts | raw counts |
| `layers["normalized"]` | log-CP10k: what the tools score | log-CP10k |
| `obs["label"]` | 29 cell types | cell types |
| `obs["orig.ident"]` | 10 patients | – |
| `uns["ccc_target"]` | `ligand, target, response` (812 rows, 166 positive) | `source, target, response` |

Sources: Wu et al. 2021 (GEO GSE176078), processed by Open Problems / Dimitrov et al. 2022;
mouse V1 (GEO GSE71585). The TNBC file is also on figshare (article 20338536), which blocks
scripted downloads.

```python
from ccc_nulls.datasets import load_openproblems
adata, truth = load_openproblems("tnbc_wu2021")   # X = log-CP10k, cell types with < 5 cells dropped
```
