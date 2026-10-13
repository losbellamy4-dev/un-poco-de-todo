#!/usr/bin/env python3
"""Run kohnert.py over a target list in parallel.  Lines of the list: q k n d.
usage: kbatch.py targets.txt secs_per_target workers log [seed_offset [extra kohnert.py args...]]
"""
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

targets = [tuple(map(int, line.split())) for line in open(sys.argv[1]) if line.strip()]
secs, workers, out = float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
seed_off = int(sys.argv[5]) if len(sys.argv) > 5 else 0
extra = sys.argv[6:]
here = os.path.dirname(os.path.abspath(__file__))


def job(t):
    q, k, n, d = t
    try:
        r = subprocess.run([sys.executable, os.path.join(here, "kohnert.py"), str(q), str(k), str(n), str(d),
                            "--secs", str(secs), "--mip-secs", "15", "--ordmin", "3", "--ordmax", "400",
                            "--seed", str(q * 100000 + k * 1000 + n + seed_off)] + extra,
                           capture_output=True, text=True, timeout=secs * 3 + 300, cwd=here)
        return t, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return t, "", "TIMEOUT"


with ThreadPoolExecutor(workers) as ex, open(out, "a") as fo:
    for fu in as_completed([ex.submit(job, t) for t in targets]):
        t, so, se = fu.result()
        fo.write(f"{t} {se.strip()} {so.strip()}\n")
        fo.flush()
        if "FOUND" in so:
            print(so.strip(), flush=True)
