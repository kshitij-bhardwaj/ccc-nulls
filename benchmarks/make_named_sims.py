"""Freeze the shared simulated datasets every tool runs on.

LR genes come from LIANA's consensus resource (real HGNC symbols, simple pairs only), so
the native tools can map them. The remaining genes are synthetic decoys. See
``ccc_nulls.simulate.simulate_named``.

Export the resource once (needs liana, in its own venv):
    python -c "import liana as li; li.resource.select_resource('consensus').to_csv(r'data/resources/liana_consensus.csv', index=False)"

Then:
    python benchmarks/make_named_sims.py                       # S0, S1, S2 x seeds 0-4
    python benchmarks/make_named_sims.py --seeds 1 --quick     # one small set for smoke tests

Writes data/sim/named_<regime>_seed<k>.h5ad plus data/sim/manifest.json. Each file has
X (log1p CP10k), obs["label"] and the truth tables in uns (read them back with
``ccc_nulls.simulate.sim_from_anndata``).
"""
import argparse
import hashlib
import json
from pathlib import Path

from ccc_nulls.simulate import load_pairs, sim_to_anndata, simulate_named

ROOT = Path(__file__).resolve().parents[1]
REGIMES = {
    "S0": dict(structured=False, planted=False, shuffle_labels=True),   # no structure, no CCC
    "S1": dict(structured=True, planted=False),                         # structure, no CCC
    "S2": dict(structured=True, planted=True),                          # structure + planted CCC
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resource", default=ROOT / "data" / "resources" / "liana_consensus.csv", type=Path)
    ap.add_argument("--out", default=ROOT / "data" / "sim", type=Path)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--n-background", type=int, default=5000)
    ap.add_argument("--n-planted", type=int, default=60)
    ap.add_argument("--plant-lfc", type=float, default=1.0)
    ap.add_argument("--quick", action="store_true", help="smaller data: 1,200 cells, 2,000 background genes")
    args = ap.parse_args()

    pairs = load_pairs(args.resource)
    kw = dict(n_background=args.n_background, n_planted=args.n_planted, plant_lfc=args.plant_lfc)
    if args.quick:
        kw.update(n_background=2000, sizes=(300, 250, 200, 200, 150, 100), n_planted=min(args.n_planted, 30))
    args.out.mkdir(parents=True, exist_ok=True)
    manifest = {"resource": str(args.resource), "resource_md5": hashlib.md5(args.resource.read_bytes()).hexdigest(),
                "n_pairs": len(pairs), "files": []}
    for regime, flags in REGIMES.items():
        for seed in range(args.seeds):
            sim = simulate_named(pairs, seed=seed, **kw, **flags)
            path = args.out / f"named_{regime}_seed{seed}.h5ad"
            sim_to_anndata(sim, regime=regime, seed=seed, **kw, **flags).write_h5ad(path)
            manifest["files"].append(dict(path=path.name, regime=regime, seed=seed, cells=len(sim.labels),
                                          genes=len(sim.genes), n_true=len(sim.truth)))
            print(f"{path.name}: {len(sim.labels)} cells, {len(sim.genes)} genes, {len(sim.truth)} planted", flush=True)
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
