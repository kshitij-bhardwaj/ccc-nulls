# AMENDMENT_9 results: does the top-1% cutoff decide the result?

Plan: [`AMENDMENT_9.md`](../../../liana_final/AMENDMENT_9.md). AMENDMENT_8's runs (same seeds, shuffles and planted pairs), scored at the top 0.5%, 1%, 2%, 5% and 10% of each ranking and by `auroc_false`, which uses no cutoff.

## Checks against AMENDMENT_8 (the 1% column)

- 897 values compared with AMENDMENT_8; 0 differ.

## False share of the top c: range over scenarios

Each cell: lowest–highest mean false share (3 replicates) over the scenarios. ✓: below liana6's in every scenario (criterion 1 at that cutoff); ✗: not.

### Mammary and Liver, both sides (4 scenarios)

| Variant | 0.5% | 1% | 2% | 5% | 10% |
|---|---:|---:|---:|---:|---:|
| liana6 | 76–86% | 80–86% | 84–89% | 84–92% | 88–94% |
| swap_null_b | 45–57% ✓ | 60–69% ✓ | 71–81% ✓ | 81–89% ✓ | 87–93% ✗ |
| swap3_null_b | 38–45% ✓ | 32–60% ✓ | 59–77% ✓ | 75–88% ✓ | 85–93% ✓ |
| swap4_null_b | 7–43% ✓ | 39–55% ✓ | 55–75% ✓ | 76–88% ✓ | 86–93% ✗ |
| swap_b2 | 60–69% ✓ | 69–78% ✓ | 75–85% ✓ | 82–92% ✗ | 87–94% ✗ |
| swap3_b2 | 51–67% ✓ | 63–78% ✓ | 70–86% ✓ | 80–92% ✗ | 87–95% ✗ |
| swap4_b2 | 51–72% ✓ | 58–80% ✓ | 67–87% ✓ | 81–93% ✗ | 87–95% ✗ |
| swap_t2 | 75–87% ✗ | 80–89% ✗ | 84–92% ✗ | 87–94% ✗ | 89–95% ✗ |
| swap3_t2 | 77–91% ✗ | 82–93% ✗ | 85–94% ✗ | 88–96% ✗ | 90–96% ✗ |
| swap4_t2 | 84–95% ✗ | 87–95% ✗ | 88–96% ✗ | 90–97% ✗ | 91–97% ✗ |

### Held-out Lung, both sides

| Variant | 0.5% | 1% | 2% | 5% | 10% |
|---|---:|---:|---:|---:|---:|
| liana6 | 79–88% | 83–88% | 90–90% | 90–93% | 91–94% |
| swap_null_b | 50–57% ✓ | 68–68% ✓ | 77–79% ✓ | 85–87% ✓ | 90–92% ✓ |
| swap3_null_b | 38–45% ✓ | 43–49% ✓ | 67–70% ✓ | 81–84% ✓ | 89–92% ✓ |
| swap4_null_b | 8–16% ✓ | 45–45% ✓ | 60–65% ✓ | 80–84% ✓ | 89–91% ✓ |
| swap_b2 | 62–69% ✓ | 74–76% ✓ | 81–83% ✓ | 87–90% ✓ | 90–94% ✓ |
| swap3_b2 | 59–63% ✓ | 70–74% ✓ | 78–82% ✓ | 86–90% ✓ | 90–94% ✓ |
| swap4_b2 | 62–67% ✓ | 71–75% ✓ | 79–82% ✓ | 86–90% ✓ | 91–94% ✓ |
| swap_t2 | 86–87% ✗ | 87–89% ✗ | 89–92% ✗ | 91–94% ✗ | 92–94% ✗ |
| swap3_t2 | 85–88% ✗ | 88–91% ✗ | 90–93% ✗ | 92–95% ✗ | 92–96% ✗ |
| swap4_t2 | 90–91% ✗ | 93–94% ✗ | 93–95% ✗ | 93–96% ✗ | 94–96% ✗ |

### Mouse brain, both sides (reported, not decided)

| Variant | 0.5% | 1% | 2% | 5% | 10% |
|---|---:|---:|---:|---:|---:|
| liana6 | 81–93% | 84–93% | 84–92% | 85–92% | 87–91% |
| swap_null_b | 64–76% ✓ | 73–82% ✓ | 78–86% ✓ | 83–87% ✓ | 88–92% ✗ |
| swap3_null_b | 31–36% ✓ | 47–53% ✓ | 67–72% ✓ | 81–84% ✓ | 87–90% ✗ |
| swap4_null_b | 12–14% ✓ | 43–45% ✓ | 64–64% ✓ | 81–83% ✓ | 88–89% ✗ |
| swap_b2 | 66–83% ✓ | 73–85% ✓ | 78–87% ✓ | 83–90% ✓ | 88–93% ✗ |
| swap3_b2 | 41–66% ✓ | 56–74% ✓ | 70–83% ✓ | 82–89% ✓ | 87–92% ✗ |
| swap4_b2 | 28–57% ✓ | 51–71% ✓ | 68–81% ✓ | 82–88% ✓ | 87–91% ✗ |
| swap_t2 | 84–93% ✗ | 87–94% ✗ | 89–94% ✗ | 88–93% ✗ | 88–92% ✗ |
| swap3_t2 | 89–94% ✗ | 90–95% ✗ | 90–95% ✗ | 88–93% ✗ | 88–92% ✗ |
| swap4_t2 | 93–95% ✗ | 93–95% ✗ | 92–94% ✗ | 91–93% ✗ | 90–92% ✗ |

## Verdicts for BPT

| Variant | Mammary and Liver | Held-out Lung |
|---|---|---|
| swap_null_b (2 tools) | holds up to 5%; fails at 10% | robust to the cutoff |
| swap3_null_b (3 tools) | robust to the cutoff | robust to the cutoff |
| swap4_null_b (4 tools) | holds up to 5%; fails at 10% | robust to the cutoff |

For MDT and LOO-t, the ✓/✗ marks above apply the same criterion (reported, not decided).

## What changes away from 1%

Every case where criterion 1 at another cutoff disagrees with its 1% result. The share of false rows among all rows (which every variant reaches at 100%) is mammary receptor 92%, mammary ligand 97%, liver receptor 93%, liver ligand 96%, lung receptor 94%, lung ligand 97%.

- Mammary and Liver, swap_null_b at 10%: fails (at 1% it passes); not below liana6 in mammary receptor (88.2% vs liana6 87.6%).
- Mammary and Liver, swap4_null_b at 10%: fails (at 1% it passes); not below liana6 in mammary receptor (88.2% vs liana6 87.6%).
- Mammary and Liver, swap_b2 at 5%: fails (at 1% it passes); not below liana6 in mammary receptor (86.5% vs liana6 84.1%).
- Mammary and Liver, swap_b2 at 10%: fails (at 1% it passes); not below liana6 in mammary receptor (89.3% vs liana6 87.6%), mammary ligand (94.3% vs liana6 93.8%).
- Mammary and Liver, swap3_b2 at 5%: fails (at 1% it passes); not below liana6 in mammary receptor (84.3% vs liana6 84.1%), mammary ligand (92.5% vs liana6 91.6%).
- Mammary and Liver, swap3_b2 at 10%: fails (at 1% it passes); not below liana6 in mammary receptor (89.5% vs liana6 87.6%), mammary ligand (95.1% vs liana6 93.8%).
- Mammary and Liver, swap4_b2 at 5%: fails (at 1% it passes); not below liana6 in mammary receptor (84.7% vs liana6 84.1%), mammary ligand (93.4% vs liana6 91.6%).
- Mammary and Liver, swap4_b2 at 10%: fails (at 1% it passes); not below liana6 in mammary receptor (89.8% vs liana6 87.6%), mammary ligand (95.4% vs liana6 93.8%).

## No cutoff: `auroc_false`

P(a row that is not false outranks a false row); higher is better. Mean over replicates, per scenario.

| Variant | mammary receptor | mammary ligand | liver receptor | liver ligand | lung receptor | lung ligand | allen brain atlas receptor | allen brain atlas ligand | above liana6 everywhere |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| liana6 | 0.498 | 0.578 | 0.519 | 0.539 | 0.514 | 0.525 | 0.496 | 0.461 |  |
| swap_null_b | 0.486 | 0.566 | 0.499 | 0.539 | 0.495 | 0.517 | 0.476 | 0.435 | no |
| swap3_null_b | 0.469 | 0.541 | 0.490 | 0.529 | 0.481 | 0.506 | 0.460 | 0.424 | no |
| swap4_null_b | 0.473 | 0.575 | 0.497 | 0.530 | 0.487 | 0.512 | 0.485 | 0.465 | no |
| swap_b2 | 0.461 | 0.546 | 0.487 | 0.504 | 0.470 | 0.480 | 0.453 | 0.416 | no |
| swap3_b2 | 0.450 | 0.514 | 0.478 | 0.480 | 0.453 | 0.465 | 0.446 | 0.420 | no |
| swap4_b2 | 0.430 | 0.467 | 0.465 | 0.450 | 0.432 | 0.441 | 0.436 | 0.412 | no |
| swap_t2 | 0.496 | 0.551 | 0.510 | 0.526 | 0.504 | 0.511 | 0.495 | 0.465 | no |
| swap3_t2 | 0.491 | 0.520 | 0.501 | 0.501 | 0.495 | 0.491 | 0.505 | 0.493 | no |
| swap4_t2 | 0.484 | 0.492 | 0.481 | 0.481 | 0.482 | 0.479 | 0.502 | 0.492 | no |

## Planted pairs (A2.2's design, f = 1.25, Mammary, 5 replicates): power in the top c, single pairs

Reported, not decided. AUROC uses no cutoff.

| Variant | 0.5% | 1% | 2% | 5% | 10% | AUROC |
|---|---:|---:|---:|---:|---:|---:|
| liana6 | 61% | 73% | 97% | 100% | 100% | 0.915 |
| swap_null_b | 92% | 99% | 100% | 100% | 100% | 0.985 |
| swap3_null_b | 96% | 99% | 100% | 100% | 100% | 0.991 |
| swap4_null_b | 97% | 99% | 100% | 100% | 100% | 0.992 |
| swap_b2 | 71% | 83% | 91% | 98% | 99% | 0.944 |
| swap3_b2 | 69% | 72% | 74% | 80% | 89% | 0.898 |
| swap4_b2 | 70% | 71% | 72% | 72% | 72% | 0.827 |
| swap_t2 | 78% | 89% | 98% | 99% | 99% | 0.952 |
| swap3_t2 | 78% | 85% | 93% | 97% | 99% | 0.937 |
| swap4_t2 | 72% | 85% | 89% | 94% | 94% | 0.910 |
