# One Kaggle session of the ccc-nulls benchmark. kaggle/kernels.py writes SLUG and COMMANDS into a copy of this file.
#
# Inputs (mounted under /kaggle/input): the code (dataset ccc-nulls-code, a `git archive` of one commit), the
# auxiliary files (dataset ccc-nulls-aux, md5-checked against benchmarks/liana_final/aux_files.md5), the five
# h5ad files (output of the session ccc-nulls-data, md5-checked when it downloaded them) and the outputs of the
# sessions this one depends on. Outputs (/kaggle/working/out): results/ (benchmarks/results/liana_final, with
# the upstream results merged in), rows/ (tables written here, for the sessions downstream) and session.json.
SLUG = "__SLUG__"
COMMANDS = __COMMANDS__
DOWNLOAD = __DOWNLOAD__          # the data session: download the h5ad files into its output instead

import glob
import json
import os
import platform
import shutil
import subprocess
import time

T0 = time.time()
OUT = "/kaggle/working/out"
CODE, DATA, PY = "/tmp/code", "/kaggle/working/data" if DOWNLOAD else "/tmp/data", "/tmp/venv/bin/python"


def sh(cmd, **kw):
    print(f"[{time.time() - T0:7.0f} s] $ {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", **kw)


def find(pattern):
    return sorted(glob.glob(f"/kaggle/input/**/{pattern}", recursive=True))


# ---------------------------------------------------------------- code and environment
tgz = find("ccc-nulls-code.tgz")
if tgz:
    os.makedirs(CODE)
    sh(f"tar xzf {tgz[0]} -C {CODE}")
else:                                   # Kaggle may have unpacked the archive
    shutil.copytree(os.path.dirname([p for p in find("pyproject.toml") if "/out/" not in p][0]), CODE)
shutil.rmtree(f"{CODE}/benchmarks/results/liana_final", ignore_errors=True)   # reproduce from scratch, not the laptop's files
commit = open(f"{CODE}/COMMIT").read().strip() if os.path.exists(f"{CODE}/COMMIT") else "unknown"
sh("pip install -q uv")
sh(f"uv venv {os.path.dirname(os.path.dirname(PY))} --python 3.12 --python-preference only-managed -q")
sh(f"uv pip install -q --python {PY} -r {CODE}/requirements/lock-liana-final.txt")
sh(f"uv pip install -q --python {PY} --no-deps -e {CODE}")
env = dict(os.environ, CCC_NULLS_DATA=DATA, PYTHONPATH=f"{CODE}:{CODE}/src", PYTHONUNBUFFERED="1",
           PATH=f"{os.path.dirname(PY)}:{os.environ['PATH']}")
os.makedirs(f"{DATA}/liana_final", exist_ok=True)

# ---------------------------------------------------------------- data
if DOWNLOAD:
    # Open Problems disabled its bucket (s3://openproblems-data, checked 2026-11-05): TNBC and Allen brain come from
    # the mirror dataset ccc-nulls-openproblems and must match the md5 in ccc_nulls.datasets; Tabula Sapiens from CELLxGENE.
    for f in find("*_log_cp10k.h5ad"):
        shutil.copy(f, DATA)
    sh(f"{PY} -m ccc_nulls.datasets tnbc_wu2021 allen_brain_atlas tabula_sapiens_mammary tabula_sapiens_liver "
       "tabula_sapiens_lung", cwd=CODE, env=env)
else:
    aux = find("ccc-nulls-aux.tgz")
    if aux:
        sh(f"tar xzf {aux[0]} -C {DATA}")
    else:
        root = os.path.dirname(os.path.dirname(find("cellphonedb_v5.0.0/cellphonedb.zip")[0]))
        shutil.copytree(root, DATA, dirs_exist_ok=True)
    sh(f"md5sum --quiet -c {CODE}/benchmarks/liana_final/aux_files.md5", cwd=DATA)
    for f in find("data/*.h5ad"):                  # the data session's output
        os.symlink(f, f"{DATA}/{os.path.basename(f)}")
    for up in sorted(glob.glob("/kaggle/input/**/out/session.json", recursive=True)):
        d = os.path.dirname(up)
        machine = json.load(open(up))["machine"]
        if os.path.isdir(f"{d}/rows"):
            shutil.copytree(f"{d}/rows", f"{DATA}/liana_final", dirs_exist_ok=True)
        sh(f"{PY} -m benchmarks.liana_final.merge_results {d}/results --machine '{machine}'", cwd=CODE, env=env)
print("machine:", platform.platform(), os.cpu_count(), "CPUs;", open("/proc/meminfo").readline().split()[1], "kB;",
      "code", commit, flush=True)

# ---------------------------------------------------------------- run
start = time.time()
os.makedirs(OUT, exist_ok=True)
status = "complete"
try:
    for c in COMMANDS:          # under kaggle/watch.py: memory every minute, the stack every 5 minutes
        sh(c.replace("python -m ", f"python {CODE}/kaggle/watch.py ", 1), cwd=CODE, env=env)
except subprocess.CalledProcessError as e:
    status = f"failed: {e.cmd}"
    raise
finally:
    if not DOWNLOAD:
        shutil.copytree(f"{CODE}/benchmarks/results/liana_final", f"{OUT}/results", dirs_exist_ok=True)
        for f in glob.glob(f"{DATA}/liana_final/**/*", recursive=True):
            if os.path.isfile(f) and not os.path.islink(f) and os.path.getmtime(f) >= start:
                dst = f"{OUT}/rows/{os.path.relpath(f, f'{DATA}/liana_final')}"
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                shutil.copy(f, dst)
        if os.path.isdir(f"{CODE}/docs/_static/figures"):
            shutil.copytree(f"{CODE}/docs/_static/figures", f"{OUT}/figures", dirs_exist_ok=True)
    owner = os.environ.get("KAGGLE_USERNAME", "kshitijbhardwaj1729")
    json.dump({"slug": SLUG, "machine": f"https://www.kaggle.com/code/{owner}/{SLUG}", "code_commit": commit,
               "commands": COMMANDS, "status": status, "seconds": round(time.time() - T0),
               "platform": platform.platform(), "cpus": os.cpu_count(),
               "mem_kb": int(open("/proc/meminfo").readline().split()[1])}, open(f"{OUT}/session.json", "w"), indent=1)
    print(f"done in {time.time() - T0:.0f} s:", status, flush=True)
