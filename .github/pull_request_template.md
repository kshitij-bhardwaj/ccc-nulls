## What and why

<!-- One or two sentences. Link the plan step (docs/PLAN.md) if relevant. -->

## Area

- [ ] `tools/<tool>/` (owner only)
- [ ] shared library `src/ccc_nulls/`
- [ ] benchmarks / results
- [ ] docs

## Checks

- [ ] `pytest` passes locally
- [ ] adapters return `ccc_nulls.schema.validate`-clean tables
- [ ] no data, large result files or personal names committed
- [ ] numbers quoted in docs come from a file in `benchmarks/results/` (path given here)
