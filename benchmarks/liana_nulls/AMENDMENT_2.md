# Amendment 2: T2 power under each tool's own call rule

Added on 1 November 2026, after `PLAN.md` and `AMENDMENT_1.md`. This is an addition, not part of the
pre-registered plan.

## What prompted it, honestly stated

A first T2 run was stopped after one cell: replicate 0, fold change 1.5. That cell's result had been
written and was then seen:

| | Null A, BH | Null B, BH |
|---|---|---|
| Power | 1.00 | 0.00 |
| FDR | 0.82 (CellChat) / 0.94 (CellPhoneDB) | 0 (no calls at all) |

The zero is structural, not a property of Null B's logic.
- With B = 1,000 the smallest p-value is 1/1001.
- BH at α = 0.05 over a family of m rows can reject only if about m / 50 rows sit at that floor.
- With 30 planted true edges, Null B (which needs both sides specific) cannot reach that.
- Null A can, because its one-sided false rows also sit at the floor.

So the pre-registered T2 criterion ("Null B's FDR upper bound ≤ 0.05") can be met by making no calls.

## What is added

1. **Tool's rule for both nulls.** In T2 each null is also scored with the tool's own unadjusted rule:
   - CellChat: p < 0.05 and prob > 0;
   - CellPhoneDB: p ≤ 0.05 among rows passing its expression gate.

   Columns `null_a_tool_rule` and `null_b_tool_rule`. Both use the repository's p convention,
   (1 + k) / (B + 1). Power, FDR and one-sided false calls are reported for both, beside the
   pre-registered BH columns.
2. **Vacuous passes are flagged.** The report marks the pre-registered T2 criterion **vacuous**
   at any effect size where Null B (BH) makes no calls on average.
3. **Extended T2 line (labelled Amendment 2).**
   - Null B's tool-rule FDR is compared with Null A's tool-rule FDR at every effect size.
   - Null B's tool-rule power is reported next to Null A's.
   - No pass/fail threshold is added after the fact. The numbers are reported as they come.

The pre-registered decision is computed exactly as written; nothing here changes it.
