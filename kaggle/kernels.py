"""The whole benchmark as Kaggle sessions (CPU, 32 GB): build, push in dependency order, collect.

    python kaggle/kernels.py code       # git archive of HEAD -> dataset ccc-nulls-code (new version)
    python kaggle/kernels.py build      # kaggle/build/<slug>/ : run.py + kernel-metadata.json per session
    python kaggle/kernels.py dispatch   # push each session once its inputs are done (at most 5 at once)
    python kaggle/kernels.py collect    # download outputs, merge into benchmarks/results/liana_final

Needs the Kaggle CLI, signed in. A session's ``needs`` are mounted as inputs; every session also mounts the code,
the auxiliary files and the h5ad files. The sessions and their links are listed in kaggle/README.md.
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

OWNER = "kshitijbhardwaj1729"
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BUILD, OUTPUTS = HERE / "build", HERE / "outputs"
MAX_ACTIVE = 5                     # Kaggle's limit on concurrent batch CPU sessions
LF, A2, A3, A4 = ("python -m benchmarks.liana_final", "python -m benchmarks.liana_final.amendment2",
                  "python -m benchmarks.liana_final.amendment3", "python -m benchmarks.liana_final.amendment4")
RUNS = {"tnbc": "tnbc_wu2021", "allen": "allen_brain_atlas", "mammary": "tabula_sapiens_mammary",
        "liver": "tabula_sapiens_liver", "hypothalamus": "moffitt2018_hypothalamus"}
LUNG = "tabula_sapiens_lung"


def S(cmds, needs=()):
    return {"commands": cmds if isinstance(cmds, list) else [cmds], "needs": list(needs)}


SESSIONS = {"data": S([], ()) | {"download": True}}
# PLAN.md + AMENDMENT_1.md
SESSIONS |= {f"run-{k}": S(f"{LF} --part run --dataset {v}", ["data"]) for k, v in RUNS.items()}
SESSIONS |= {f"calibration-{k}": S(f"{LF} --part calibration --dataset {RUNS[k]}", ["data"]) for k in ("mammary", "liver")}
SESSIONS["power"] = S(f"{LF} --part power", ["data"])
SESSIONS["biology"] = S(f"{LF} --part biology", ["data"] + [f"run-{k}" for k in RUNS])
# AMENDMENT_2.md
for scen in ("receptor", "ligand"):
    for rep in range(3):
        SESSIONS[f"a2-lung-{scen}-{rep}"] = S(f"{A2} --part confirm --only '^{LUNG}\\|{scen}\\|{rep}$'", ["data"])
SESSIONS["a2-lung-power"] = S(f"{A2} --part confirm --only '^{LUNG}\\|f2\\|'", ["data"])
for scheme in ("global", "donor"):
    for rep in range(3):
        SESSIONS[f"a2-lung-donor-{scheme}-{rep}"] = S(f"{A2} --part donor --only '^{LUNG}\\|{scheme}\\|{rep}$'", ["data"])
SESSIONS["a2-hard-power"] = S(f"{A2} --part hard_power", ["data"])
SESSIONS |= {f"a2-recal-{k}": S(f"{A2} --part recal --only '^{RUNS[k]}\\|'", ["data"]) for k in ("mammary", "liver")}
SESSIONS["a2-mouse"] = S(f"{A2} --part mouse", ["data"])
SESSIONS["a2-timing"] = S(f"{A2} --part timing", ["data"])
SESSIONS |= {f"a2-donor-{k}": S(f"{A2} --part donor --only '^{RUNS[k]}\\|'", ["data"]) for k in ("tnbc", "mammary", "liver")}
SESSIONS["a2-donor-real"] = S(f"{A2} --part donor --only '^({RUNS['tnbc']}|{RUNS['mammary']}|{RUNS['liver']})$'",
                              ["data", "run-tnbc", "run-mammary", "run-liver"])
SESSIONS["a2-sensitivity-seed-expr"] = S(f"{A2} --part sensitivity --only '^seed|@expr_prop='", ["data", "run-mammary"])
SESSIONS["a2-sensitivity-cells-labels"] = S(f"{A2} --part sensitivity --only '@(cells|labels)='", ["data", "run-mammary"])
SESSIONS["a2-sensitivity-perms-hypothalamus"] = S(f"{A2} --part sensitivity --only '^real$|^{RUNS['hypothalamus']}@n_perms'",
                                                  ["data", "run-mammary", "run-hypothalamus"])
SESSIONS["a2-sensitivity-perms-mammary"] = S(f"{A2} --part sensitivity --only '^real@|^{RUNS['mammary']}@n_perms'",
                                             ["data", "run-mammary", "run-hypothalamus"])
SESSIONS["a2-cases"] = S(f"{A2} --part cases", ["data"] + [f"run-{k}" for k in RUNS])
SESSIONS["a2-biology"] = S(f"{A2} --part biology", ["data", "biology", "a2-donor-real"] + [f"run-{k}" for k in RUNS])
# AMENDMENT_3.md: T2 (PR #11) in the same runs
SESSIONS |= {f"a3-real-{k}": S(f"{A3} --part real --only '^{v}$'", ["data", f"run-{k}"]) for k, v in RUNS.items()}
SESSIONS |= {f"a3-calibration-{k}": S(f"{A3} --part calibration --only '^{RUNS[k]}\\|'", ["data"]) for k in ("mammary", "liver")}
SESSIONS["a3-power"] = S(f"{A3} --part power", ["data"])
SESSIONS["a3-hard-power"] = S(f"{A3} --part hard_power", ["data"])
SESSIONS["a3-biology"] = S(f"{A3} --part biology", ["data", "biology"] + [f"a3-real-{k}" for k in RUNS])
# AMENDMENT_4.md: T2 on the rest of AMENDMENT_2
for scen in ("receptor", "ligand"):
    for rep in range(3):
        SESSIONS[f"a4-lung-{scen}-{rep}"] = S(f"{A4} --part lung --only '^{LUNG}\\|{scen}\\|{rep}$'", ["data"])
SESSIONS["a4-lung-power"] = S(f"{A4} --part lung --only '^{LUNG}\\|f2\\|'", ["data"])
SESSIONS |= {f"a4-donor-lung-{rep}": S(f"{A4} --part donor --only '^{LUNG}\\|global\\|{rep}$'", ["data"]) for rep in range(3)}
SESSIONS |= {f"a4-donor-{k}": S(f"{A4} --part donor --only '^{RUNS[k]}\\|'", ["data"]) for k in ("tnbc", "mammary", "liver")}
SESSIONS["a4-donor-real"] = S(f"{A4} --part donor --only '^({RUNS['tnbc']}|{RUNS['mammary']}|{RUNS['liver']})$'",
                              ["data", "a3-real-tnbc", "a3-real-mammary", "a3-real-liver"])
SESSIONS["a4-sensitivity-seed-expr"] = S(f"{A4} --part sensitivity --only '^seed|@expr_prop='", ["data", "a3-real-mammary"])
SESSIONS["a4-sensitivity-cells-labels"] = S(f"{A4} --part sensitivity --only '@(cells|labels)='", ["data", "a3-real-mammary"])
SESSIONS["a4-sensitivity-perms-hypothalamus"] = S(f"{A4} --part sensitivity --only '^real$|^{RUNS['hypothalamus']}@n_perms'",
                                                  ["data", "a3-real-hypothalamus"])
SESSIONS["a4-sensitivity-perms-mammary"] = S(f"{A4} --part sensitivity --only '^real@|^{RUNS['mammary']}@n_perms'",
                                             ["data", "a3-real-mammary"])
SESSIONS["a4-mouse"] = S(f"{A4} --part mouse", ["data"])
SESSIONS["a4-cases"] = S(f"{A4} --part cases", ["data"] + [f"run-{k}" for k in RUNS] + [f"a3-real-{k}" for k in RUNS])
SESSIONS["a4-timing"] = S(f"{A4} --part timing", ["data", "run-mammary", "run-hypothalamus"])
# everything merged; reports and figures
SESSIONS["report"] = S([f"{LF} --part report", f"{A2} --part report", f"{A3} --part report", f"{A4} --part report",
                        "python -m benchmarks.liana_final.figures"],
                       [k for k in SESSIONS if k != "data"])


def slug(k: str) -> str:
    return f"ccc-nulls-{k}"


def kaggle(*args, check=True) -> str:
    r = subprocess.run(["kaggle", *args], capture_output=True, text=True)
    if check and r.returncode:
        raise SystemExit(r.stdout + r.stderr)
    return (r.stdout + r.stderr).strip()


def code():
    d = BUILD / "_code"
    d.mkdir(parents=True, exist_ok=True)
    sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()
    if subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, capture_output=True, text=True).stdout:
        raise SystemExit("commit first: the sessions run a committed version of the code")
    subprocess.run(f"git archive HEAD --prefix=ccc-nulls/ -o {d}/code.tar && mkdir -p {d}/x && tar xf {d}/code.tar -C {d}/x "
                   f"&& echo {sha} > {d}/x/ccc-nulls/COMMIT && tar czf {d}/ccc-nulls-code.tgz -C {d}/x/ccc-nulls . "
                   f"&& rm -rf {d}/x {d}/code.tar", shell=True, check=True)
    (d / "dataset-metadata.json").write_text(json.dumps({"title": "ccc-nulls-code", "id": f"{OWNER}/ccc-nulls-code",
                                                         "licenses": [{"name": "other"}]}))
    print(kaggle("datasets", "version", "-p", str(d), "-m", f"code at {sha}", check=False))


def build():
    template = (HERE / "runner.py").read_text()
    for k, s in SESSIONS.items():
        d = BUILD / slug(k)
        d.mkdir(parents=True, exist_ok=True)
        (d / "run.py").write_text(template.replace("__SLUG__", slug(k)).replace("__COMMANDS__", json.dumps(s["commands"]))
                                  .replace("__DOWNLOAD__", str(bool(s.get("download")))))
        (d / "kernel-metadata.json").write_text(json.dumps({
            "id": f"{OWNER}/{slug(k)}", "title": slug(k), "code_file": "run.py", "language": "python",
            "kernel_type": "script", "is_private": True, "enable_gpu": False, "enable_tpu": False, "enable_internet": True,
            "dataset_sources": [f"{OWNER}/ccc-nulls-code"] + ([f"{OWNER}/ccc-nulls-openproblems"] if s.get("download")
                                                             else [f"{OWNER}/ccc-nulls-aux"]),
            "kernel_sources": [f"{OWNER}/{slug(n)}" for n in s["needs"]], "competition_sources": [], "model_sources": []},
            indent=1))
    print(f"built {len(SESSIONS)} sessions in {BUILD}")


def status(k: str) -> str:
    s = kaggle("kernels", "status", f"{OWNER}/{slug(k)}", check=False)
    for w in ("COMPLETE", "ERROR", "CANCEL", "RUNNING", "QUEUED"):
        if w in s:
            return w.lower()
    return "new"


def dispatch(only=None):
    """Push sessions whose inputs are complete; a session already complete is not pushed again."""
    todo = [k for k in SESSIONS if only is None or k in only]
    state = {k: status(k) for k in SESSIONS}
    state |= {k: "new" for k in todo if state[k] in ("error", "cancel")}       # a failed session runs again
    while True:
        for k in todo:
            if state[k] in ("running", "queued"):
                state[k] = status(k)
                if state[k] == "complete":
                    print(time.strftime("%H:%M"), "complete", k, flush=True)
                elif state[k] in ("error", "cancel"):
                    print(time.strftime("%H:%M"), "FAILED", k, flush=True)
                    return 1
        active = sum(state[k] in ("running", "queued") for k in SESSIONS)
        ready = [k for k in todo if state[k] == "new" and all(state[n] == "complete" for n in SESSIONS[k]["needs"])]
        for k in ready[:max(0, MAX_ACTIVE - active)]:
            r = kaggle("kernels", "push", "-p", str(BUILD / slug(k)), check=False)
            if "successfully" in r:
                state[k] = "running"
                print(time.strftime("%H:%M"), "pushed", k, flush=True)
            else:
                print(time.strftime("%H:%M"), "push refused", k, r.splitlines()[-1], flush=True)
                break
        if all(state[k] == "complete" for k in todo):
            print(time.strftime("%H:%M"), "ALL COMPLETE", flush=True)
            return 0
        time.sleep(120)


def collect():
    """Download every session's output and merge it into benchmarks/results/liana_final (each run keeps its link)."""
    sys.path.insert(0, str(ROOT))
    from benchmarks.liana_final.merge_results import merge
    for k in SESSIONS:
        if SESSIONS[k].get("download"):
            continue
        d = OUTPUTS / slug(k)
        if not (d / "out" / "session.json").exists():
            d.mkdir(parents=True, exist_ok=True)
            kaggle("kernels", "output", f"{OWNER}/{slug(k)}", "-p", str(d), check=False)
        meta = json.loads((d / "out" / "session.json").read_text())
        merge(d / "out" / "results", meta["machine"])


def readme():
    """kaggle/README.md: every session, its link, inputs and command (with run time once collected)."""
    url = lambda k: f"https://www.kaggle.com/code/{OWNER}/{slug(k)}"
    L = ["# The benchmark on Kaggle", "",
         "Every result in `benchmarks/results/liana_final/` was produced by one of the Kaggle sessions below (CPU, 4 cores, "
         "32 GB), from one commit of this repository, in the locked environment (`requirements/lock-liana-final.txt`). "
         "Each result file, and each run inside AMENDMENT_2's files, records the link of the session that produced it "
         "(`machine`). To rerun: `python kaggle/kernels.py code && python kaggle/kernels.py build && "
         "python kaggle/kernels.py dispatch && python kaggle/kernels.py collect`. To run elsewhere: [`RUN.md`](../RUN.md).", "",
         "Inputs:", "",
         f"- code: [ccc-nulls-code](https://www.kaggle.com/datasets/{OWNER}/ccc-nulls-code) (`git archive` of the commit in `COMMIT`)",
         f"- auxiliary files: [ccc-nulls-aux](https://www.kaggle.com/datasets/{OWNER}/ccc-nulls-aux) "
         "(md5: `benchmarks/liana_final/aux_files.md5`)",
         f"- h5ad files: the output of [ccc-nulls-data]({url('data')}): Tabula Sapiens from CELLxGENE, TNBC and Allen brain "
         f"from the [Open Problems copies](https://www.kaggle.com/datasets/{OWNER}/ccc-nulls-openproblems) (the Open Problems "
         "bucket is closed), all md5-checked", "",
         "The sessions and data sets are private to the project team; to be added as a collaborator, ask a maintainer.", "",
         "| Session | Needs | Command | Run time |", "|---|---|---|---|"]
    for k, s in SESSIONS.items():
        f = OUTPUTS / slug(k) / "out" / "session.json"
        secs = json.loads(f.read_text())["seconds"] if f.exists() else None
        cmd = "download the five h5ad files" if s.get("download") else "<br>".join(f"`{c}`" for c in s["commands"])
        needs = ", ".join(n for n in s["needs"] if n != "data") if k != "report" else "all"
        L.append(f"| [{k}]({url(k)}) | {needs or '–'} | {cmd.replace('|', '\\|')} | "
                 f"{'' if secs is None else f'{secs / 3600:.1f} h'} |")
    (HERE / "README.md").write_text("\n".join(L) + "\n")
    print("wrote", HERE / "README.md")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    {"code": code, "build": build, "dispatch": lambda: sys.exit(dispatch(sys.argv[2:] or None)), "collect": collect,
     "readme": readme}.get(cmd, lambda: sys.exit(__doc__))()
