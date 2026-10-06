"""Quickstart on simulated data with known truth.

python examples/quickstart.py
"""
import ccc_nulls as cn
from ccc_nulls.simulate import candidates_from_pairs, label_truth, simulate

# Structured data: cell identity nudges every gene, plus 60 planted LR events (LFC 1.5).
sim = simulate(seed=0, plant_lfc=1.5)
cands = candidates_from_pairs(sim.X, sim.genes, sim.labels, sim.pairs)   # common schema

res = cn.significance_layer(
    sim.X, sim.genes, sim.labels, cands,
    score="mean",                                   # CellPhoneDB statistic
    B=2000,                                         # Null A permutations
    resource_genes=set(sim.pairs.values.ravel()),   # keep LR genes out of the decoy pool
)
res["truth"] = label_truth(res, sim)

print(f"{len(res):,} candidate tests, {int((res.truth == 'true').sum())} planted")
print(cn.overlap_summary(res).to_string(), "\n")
for name, q in (("Null A", "q_rand"), ("Null B", "q_strat")):
    called = res[q] < 0.05
    tp = (called & (res.truth == "true")).sum()
    fp = (called & (res.truth == "null")).sum()
    print(f"{name}: {called.sum():5d} calls | FDR {fp / max(tp + fp, 1):.1%} | power {tp / (res.truth == 'true').sum():.0%}")
