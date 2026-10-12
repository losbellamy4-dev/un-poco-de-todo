#!/usr/bin/env python3
"""Run qcyclic.py over many [[n,k]] targets in parallel.

Targets: n in [nmin, nmax], k in [kmin, kmax], lb < ub, k reachable by the
cyclic family; ordered by lower-bound plateau length (long plateaus = bounds
obtained by padding, likely weak).
usage: qbatch.py --nmin 30 --nmax 64 --kmin 1 --kmax 20 --secs 90 --workers 4 --out runs/qb.log
"""
import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

from qcyclic import CyclicFamily
from qsearch import load_qbounds


def targets(nmin, nmax, kmin, kmax, oddonly):
    B = load_qbounds()
    fams, out = {}, []
    for (n, k), (lo, hi) in B.items():
        if not (nmin <= n <= nmax and kmin <= k <= kmax) or lo >= hi:
            continue
        if oddonly and n % 2 == 0:
            continue
        L = 0
        while (n - L - 1, k) in B and B[(n - L - 1, k)][0] == lo:
            L += 1
        if n not in fams:
            fams[n] = CyclicFamily(n)
        if k in fams[n].reach:
            out.append((L, n, k))
    out.sort(reverse=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmin", type=int, default=30)
    ap.add_argument("--nmax", type=int, default=64)
    ap.add_argument("--kmin", type=int, default=0)
    ap.add_argument("--kmax", type=int, default=20)
    ap.add_argument("--secs", type=float, default=90)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--odd", action="store_true")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    T = targets(args.nmin, args.nmax, args.kmin, args.kmax, args.odd)
    print(f"{len(T)} targets", file=sys.stderr, flush=True)
    here = os.path.dirname(os.path.abspath(__file__))

    def job(t):
        L, n, k = t
        try:
            r = subprocess.run([sys.executable, os.path.join(here, "qcyclic.py"), str(n), str(k), "--secs",
                                str(args.secs), "--seed", str(n * 1000 + k)], capture_output=True, text=True,
                               timeout=3 * args.secs + 120)
            return t, r.stdout, r.stderr
        except subprocess.TimeoutExpired as e:
            return t, (e.stdout or b"").decode() if isinstance(e.stdout, bytes) else (e.stdout or ""), "TIMEOUT"

    with ThreadPoolExecutor(args.workers) as ex, open(args.out, "a") as fo:
        for fu in as_completed([ex.submit(job, t) for t in T]):
            t, out, err = fu.result()
            fo.write(f"{t} {err.strip()}\n")
            fo.flush()
            for line in out.splitlines():
                if "CANDIDATE" in line:
                    print(line, flush=True)


if __name__ == "__main__":
    main()
