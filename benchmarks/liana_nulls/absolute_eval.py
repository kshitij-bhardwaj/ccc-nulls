"""Biology without reference to Null A (NOT pre-registered; an addition).

The plan's endpoint compares, among Null A's calls, the calls another null keeps with those it drops. Here every null,
Null A included, is compared with ALL tested rows (the tool's BH family, source != target):

    share       fraction of the null's calls with evidence
    base        fraction of all tested rows with evidence (chance level)
    enrichment  share / base
    odds_ratio  odds of evidence, called vs not called, over all tested rows; 95% CI from 2,000 bootstrap resamples
                of the evidence units (as the plan's endpoints): (source, target) pairs for spatial truths (Xenium:
                its broad class pairs), (footprint, receiver) for downstream truths.

Evidence: spatial (TNBC Visium, TNBC Xenium proxy, Allen brain, hypothalamus MERFISH) and downstream response
(CytoSig, PROGENy, NicheNet; downstream.truth_rows on every tested row).

    python -m benchmarks.liana_nulls.absolute_eval --tool natmi
Results: benchmarks/results/liana_nulls/addendum_absolute_<tool>.json
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from . import downstream
from .common import AMENDMENT_DATASETS, CORE_DATASETS, OUT, log, save
from .nichenet_eval import profiles_for, rows_path, spatial_truths

NULLS = ("null_a", "null_b", "b2", "loo_t")


def odds(k):
    a, b, c, d = k + (0.5 if (k == 0).any() else 0)
    return (a / b) / (c / d)


def compare(evidence: np.ndarray, unit: np.ndarray, called: np.ndarray, n_boot: int = 2000, seed: int = 0) -> dict:
    """Called vs not called over all rows; bootstrap over units."""
    ev, called = evidence.astype(bool), called.astype(bool)
    d = pd.DataFrame({"u": unit, "a": called & ev, "b": called & ~ev, "c": ~called & ev, "e": ~called & ~ev})
    K = d.groupby("u")[["a", "b", "c", "e"]].sum().to_numpy(float)
    tot = K.sum(0)
    out = dict(calls=int(tot[0] + tot[1]), share=float(tot[0] / max(tot[0] + tot[1], 1)), base=float(ev.mean()),
               n_rows=int(len(ev)), n_units=int(len(K)))
    out["enrichment"] = out["share"] / out["base"] if out["base"] > 0 and out["calls"] else None
    if out["calls"] == 0 or tot[2] + tot[3] == 0:
        return out | dict(odds_ratio=None, ci95=None)
    rng = np.random.default_rng(seed)
    boot = [odds(K[rng.integers(0, len(K), len(K))].sum(0)) for _ in range(n_boot)]
    return out | dict(odds_ratio=float(odds(tot)), ci95=[float(x) for x in np.percentile(boot, [2.5, 97.5])])


def main(argv=None):
    from .amendment4 import loo_t
    from .data import SPECIES
    from .xenium import TNBC_CLASS

    ap = argparse.ArgumentParser(prog="python -m benchmarks.liana_nulls.absolute_eval")
    ap.add_argument("--tool", default="natmi")
    ap.add_argument("--dataset", nargs="+", default=CORE_DATASETS + AMENDMENT_DATASETS)
    args = ap.parse_args(argv)
    name_out = f"addendum_absolute_{args.tool}"
    out = {}
    for name in args.dataset:
        rows = pd.read_pickle(rows_path(args.tool, name))
        for c in ("source", "target"):
            rows[c] = rows[c].astype(str)
        if "call_loo_t" not in rows:
            rows["p_loo_t"], rows["call_loo_t"] = loo_t(rows)
        U = rows[(rows.source != rows.target) & rows.family.to_numpy(bool)].reset_index(drop=True)
        res = {"tested_rows": int(len(U))}

        for tname, truth in spatial_truths(name).items():
            if tname == "xenium_proxy":                   # class-level truth, class-pair units
                cls = pd.read_csv(OUT / "truth_xenium_breast_classes.csv").set_index(["source_class", "target_class"]).colocalised
                s, t = U.source.map(TNBC_CLASS), U.target.map(TNBC_CLASS)
                ev = np.array([cls.get(k, np.nan) if isinstance(k[0], str) and isinstance(k[1], str) else np.nan
                               for k in zip(s, t)], float)
                unit = (s.astype(str) + "|" + t.astype(str)).to_numpy()
            else:
                tt = truth.set_index(["source", "target"]).colocalised
                ev = np.array([tt.get(k, np.nan) for k in zip(U.source, U.target)], float)
                unit = (U.source + "|" + U.target).to_numpy()
            ok = np.isfinite(ev)
            res[f"spatial_{tname}"] = {n: compare(ev[ok] > 0, unit[ok], U[f"call_{n}"].to_numpy(bool)[ok]) for n in NULLS}
            log(f"{name} spatial {tname}: " + "; ".join(f"{n} {v['enrichment'] and round(v['enrichment'], 2)}x"
                                                        for n, v in res[f'spatial_{tname}'].items()))

        prof, adata, _ = profiles_for(name)
        del adata
        rh = downstream.human_rows(U.assign(call_null_a=True), SPECIES[name])   # truth_rows keeps call_null_a rows
        for c in NULLS:
            rh[f"called_{c}"] = U[f"call_{c}"].to_numpy(bool)
        for resource in ("cytosig", "progeny", "nichenet"):
            r = downstream.truth_rows(rh, prof, resource)
            if not len(r):
                continue
            res[resource] = {n: compare(r.responds.to_numpy(bool), r.unit.to_numpy(), r[f"called_{n}"].to_numpy(bool)) for n in NULLS}
            log(f"{name} {resource}: " + "; ".join(f"{n} {v['enrichment'] and round(v['enrichment'], 2)}x"
                                                  for n, v in res[resource].items()))
        out[name] = res
        save(name_out, {"addendum": "biology without reference to Null A: each null's calls vs all tested rows",
                        "tool": args.tool, "results": out})


if __name__ == "__main__":
    main()
