"""python kaggle/watch.py <module> [args]: run ``python -m <module>`` printing memory every minute (resident,
peak, available), and the stack on a fatal signal. (Periodic stack dumps, ``faulthandler.dump_traceback_later``,
read other threads without the GIL and crashed a session with SIGSEGV; they are not used.)"""
import faulthandler
import runpy
import sys
import threading
import time


def mem() -> str:
    f = {l.split(":")[0]: int(l.split()[1]) for l in open("/proc/self/status") if l.startswith(("VmRSS", "VmHWM"))}
    avail = next(int(l.split()[1]) for l in open("/proc/meminfo") if l.startswith("MemAvailable"))
    return f"rss {f['VmRSS'] / 2**20:.1f} GB, peak {f['VmHWM'] / 2**20:.1f} GB, available {avail / 2**20:.1f} GB"


def loop():
    while True:
        time.sleep(60)
        print(time.strftime("%H:%M:%S"), "memory:", mem(), flush=True)


threading.Thread(target=loop, daemon=True).start()
faulthandler.enable()
mod, sys.argv = sys.argv[1], [sys.argv[1]] + sys.argv[2:]
runpy.run_module(mod, run_name="__main__", alter_sys=True)
