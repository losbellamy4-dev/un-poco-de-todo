#!/usr/bin/env python3
"""Multi-generator Hermitian self-orthogonal QC codes over GF(4), n <= 128.

If <g> (length l) is Hermitian self-orthogonal, every code inside g*R^m is,
so t generators (g f^(j)_1, ..., g f^(j)_m) give C = [m*l, t*kappa0] with
kappa0 = l - deg g, hence [[m*l, m*l - 2*t*kappa0, d]],  d = wt(C^perp_h \\ C).
usage: qqc2.py N K [--secs S] [--seed s]
"""
import argparse
import json
import random
import sys
import time

import numpy as np

from gf import ROOT
from qqc import qc_generator, self_orth_divisors, to_symplectic
from qsearch import load_qbounds
from quantum import normalizer, qdist128, rank


def options(N, K):
    out = []
    if (N - K) % 2:
        return out
    kap = (N - K) // 2
    for m in range(2, 17):
        if N % m:
            continue
        l = N // m
        if l % 2 == 0 or l < 3:
            continue
        try:
            divs = self_orth_divisors(l)
        except ValueError:
            continue
        for deg, gs in divs.items():
            k0 = l - deg
            if k0 and kap % k0 == 0 and 1 <= kap // k0 <= m:
                out.append((m, l, kap // k0, gs))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int)
    ap.add_argument("K", type=int)
    ap.add_argument("--secs", type=float, default=120)
    ap.add_argument("--seed", type=int, default=int(time.time()))
    args = ap.parse_args()
    N, K = args.N, args.K
    lo, hi = load_qbounds()[(N, K)]
    T = lo + 1
    opts = options(N, K)
    if not opts:
        print(f"[[{N},{K}]]: no option", file=sys.stderr)
        return
    rng = random.Random(args.seed)
    t0 = time.time()
    tried = 0
    while time.time() - t0 < args.secs:
        m, l, t, gs = rng.choice(opts)
        g = rng.choice(gs)
        G = np.vstack([qc_generator(g, [[rng.randrange(4) for _ in range(l)] for _ in range(m)], l)
                       for _ in range(t)])
        S = to_symplectic(G)
        if rank(S, N) != N - K:
            continue
        tried += 1
        Nn = normalizer(S, N)
        d = qdist128(Nn, N, target=T, threads=1, exclude=S)
        if d >= T:
            d = qdist128(Nn, N, threads=1, exclude=S)
            rec = dict(n=N, k=K, d=d, table=[lo, hi], family="herm-QC-multi", m=m, l=l, t=t,
                       stab=[list(v) for v in S], time=time.time())
            with open(ROOT / "results" / "qcandidates.jsonl", "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(f"*** QUANTUM CANDIDATE [[{N},{K},{d}]] (table {lo}-{hi}) m={m} l={l} t={t}", flush=True)
    print(f"[[{N},{K}]] herm-QC-multi: tried {tried}, {time.time()-t0:.0f}s, opts {[(m, l, t) for m, l, t, _ in opts]}",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
