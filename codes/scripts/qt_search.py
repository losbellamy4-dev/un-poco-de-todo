#!/usr/bin/env python3
"""Random search of 1-generator quasi-twisted codes against the codetables bounds.

For a field q, enumerates targets (l, m, k) with n = m*l a gap in the table and
cheap verification, then spends a time budget per target sampling random codes.
Each code is tested first against d >= lower (matching the record) and, if it
passes, against d >= lower+1 (a new record).  Hits are appended to
codes/results/candidates.jsonl for independent verification.

usage: qt_search.py q --worker i --workers W [--secs S] [--maxcost C] [--nmin N --nmax N]
"""
import argparse
import json
import random
import sys
import time

import numpy as np

from gf import ROOT, min_distance
from qt import QTFamily, bz_cost, load_bounds


def targets(q, B, maxcost, nmin, nmax, kmin, kmax, mmax):
    out = []
    for (n, k), (lo, hi) in B.items():
        if lo >= hi or not (nmin <= n <= nmax) or not (kmin <= k <= kmax):
            continue
        c = bz_cost(q, n, k, lo + 1)
        if c > maxcost:
            continue
        for m in range(2, mmax + 1):
            if n % m:
                continue
            l = n // m
            if l < k or l % GF_CHAR[q] == 0:
                continue
            out.append((c, n, k, l, m, lo, hi))
    out.sort()
    return out


GF_CHAR = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2, 9: 3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int)
    ap.add_argument("--worker", type=int, default=0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--secs", type=float, default=20)
    ap.add_argument("--maxcost", type=float, default=3e8)
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=999)
    ap.add_argument("--kmin", type=int, default=3)
    ap.add_argument("--kmax", type=int, default=999)
    ap.add_argument("--mmax", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--stats", default=None)
    args = ap.parse_args()
    q = args.q
    B = load_bounds(q)
    T = targets(q, B, args.maxcost, args.nmin, args.nmax, args.kmin, args.kmax, args.mmax)
    T = T[args.worker::args.workers]
    rng = random.Random(args.seed * 1000 + args.worker)
    fam_cache = {}
    cand_path = ROOT / "results" / "candidates.jsonl"
    stats_f = open(args.stats, "a") if args.stats else None
    print(f"q={q} worker {args.worker}: {len(T)} targets", file=sys.stderr, flush=True)
    for cost, n, k, l, m, lo, hi in T:
        lams = list(range(1, q))
        t0 = time.time()
        tried = matched = 0
        best_seen = 0
        while time.time() - t0 < args.secs:
            lam = rng.choice(lams)
            key = (l, lam)
            if key not in fam_cache:
                fam_cache[key] = QTFamily(q, l, lam)
            fam = fam_cache[key]
            if (l - k) not in fam.reach:
                lams.remove(lam)
                if not lams:
                    break
                continue
            g = fam.random_divisor(l - k, rng)
            fs = [[1]] + [[rng.randrange(q) for _ in range(k)] for _ in range(m - 1)]
            G = fam.generator(g, fs, k)
            tried += 1
            d1 = min_distance(q, G, target=lo, threads=1)
            if d1 < lo:
                best_seen = max(best_seen, d1)
                continue
            matched += 1
            d2 = min_distance(q, G, target=lo + 1, threads=1)
            if d2 >= lo + 1:
                rec = dict(q=q, n=n, k=k, d_lower_bound=d2, table=[lo, hi], kind="QT",
                           l=l, m=m, lam=lam, g=g, fs=fs, time=time.time())
                with open(cand_path, "a") as f:
                    f.write(json.dumps(rec) + "\n")
                print(f"*** CANDIDATE q={q} [{n},{k}] d>={d2} (table {lo}-{hi})", flush=True)
        line = dict(q=q, n=n, k=k, l=l, m=m, lo=lo, hi=hi, cost=cost, tried=tried,
                    matched=matched, best_below=best_seen)
        if stats_f:
            stats_f.write(json.dumps(line) + "\n")
            stats_f.flush()
        print(f"q={q} [{n},{k}] l={l} m={m} table {lo}-{hi}: tried {tried}, matched {matched}",
              file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
