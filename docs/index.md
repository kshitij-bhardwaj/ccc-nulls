# ccc-nulls

**Is the sender specific *and* the receiver specific?**

Permutation-based cell–cell communication tools (CellPhoneDB, CellChat, and LIANA+'s consensus built on them)
call an interaction "specific" as soon as **either** the ligand side or the receptor side is specific to its
cell type. ccc-nulls replaces their p-value with **Null B**, a composite null under which **both** sides must be
specific, computed from LIANA+'s own permutations, and slots it into LIANA+'s consensus ranking.

```python
import ccc_nulls as cn
res = cn.rank_aggregate(adata, groupby="cell_type")      # LIANA+'s rank_aggregate, with Null B
```

The result has LIANA+'s columns (`magnitude_rank`, `specificity_rank`, each method's scores), so LIANA+'s plots and
downstream tools work unchanged.

- **[Quickstart](notebooks/quickstart.ipynb)**: run it on bundled data in a minute.
- **[Method](method.md)**: the intersection–union test, and how it sits inside LIANA+.
- **[Benchmark](benchmark.md)**: calibration, power, biology, held-out confirmation, robustness.
- **[Reproducibility](reproducibility.md)**: pre-registration, pinned environment, one command.

```{toctree}
:hidden:
:maxdepth: 2

installation
notebooks/quickstart
method
benchmark
reproducibility
api
earlier_work
PLAN
```
