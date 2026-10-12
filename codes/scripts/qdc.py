#!/usr/bin/env python3
"""Quantum codes from 'unitary double circulant' GF(4) codes.

R = GF(4)[x]/(x^l - 1), l odd.  u = b / conj(b)^* (b invertible) satisfies
u * conj(u)^* = 1, hence for any divisor g of x^l - 1 the code
    C = {(c, c*u) : c in <g>}         (dimension l - deg g)
is Hermitian self-orthogonal, giving [[2l, 2*deg g, d]] with
d = wt(C^perp_h \\ C).  g = 1 gives Hermitian self-dual double circulant codes.
(C^perp_h contains the self-dual code {(y, y*u)}, so d is capped by its distance.)
usage: qdc.py N K [--secs S] [--seed s]
"""
import argparse
import json
import random
import sys
import time

import numpy as np

from gf import GF, ROOT
from poly import PolyRing
from qqc import conj_rev, mulmod, to_symplectic
from qsearch import load_qbounds
from quantum import normalizer, qdist, qdist128, rank

F4 = GF(4)
R4 = PolyRing(4)


def inverse_mod(b, l):
    """Inverse of b in GF(4)[x]/(x^l - 1) or None."""
    m = R4.x_n_minus(l, 1)
    r0, r1 = m, R4.trim(list(b))
    s0, s1 = [], [1]
    while r1:
        qt, rem = R4.divmod(r0, r1)
        r0, r1 = r1, rem
        s0, s1 = s1, R4.sub(s0, R4.mul(qt, s1))
    if len(r0) != 1:
        return None
    out = R4.scale(int(F4.inv[r0[0]]), s0)
    return R4.mod(out, m) if out else None


def divisors_by_degree(l):
    facs = R4.factor_squarefree(R4.x_n_minus(l, 1))
    out = {}
    nf = min(len(facs), 16)
    for mask in range(1 << nf):
        g = [1]
        for i in range(nf):
            if mask >> i & 1:
                g = R4.mul(g, facs[i])
        out.setdefault(len(g) - 1, []).append(g)
    return out


def build(l, g, u):
    k = l - (len(g) - 1)
    rows = []
    v = list(g) + [0] * (l - len(g))
    for _ in range(k):
        cu = mulmod(R4.trim(list(v)), u, l)
        rows.append(list(v) + (cu + [0] * (l - len(cu))))
        v = [v[-1]] + v[:-1]
    return np.array(rows, np.uint8)


def random_unitary(l, rng):
    while True:
        b = R4.trim([rng.randrange(4) for _ in range(l)])
        inv = inverse_mod(conj_rev(b, l), l)
        if inv is not None:
            return mulmod(b, inv, l)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int)
    ap.add_argument("K", type=int)
    ap.add_argument("--secs", type=float, default=120)
    ap.add_argument("--seed", type=int, default=int(time.time()))
    args = ap.parse_args()
    N, K = args.N, args.K
    if N % 2 or K % 2 or (N // 2) % 2 == 0:
        print(f"[[{N},{K}]]: needs N = 2l with l odd and K even", file=sys.stderr)
        return
    l = N // 2
    divs = divisors_by_degree(l).get(K // 2, [])
    if not divs:
        print(f"[[{N},{K}]]: no divisor of degree {K // 2}", file=sys.stderr)
        return
    lo, hi = load_qbounds()[(N, K)]
    T = lo + 1
    dist = qdist if N <= 64 else qdist128
    rng = random.Random(args.seed)
    t0 = time.time()
    tried = 0
    while time.time() - t0 < args.secs:
        u = random_unitary(l, rng)
        g = rng.choice(divs)
        S = to_symplectic(build(l, g, u))
        if rank(S, N) != N - K:
            continue
        tried += 1
        if K == 0:
            d = dist(S, N, target=T, threads=1)
        else:
            d = dist(normalizer(S, N), N, target=T, threads=1, exclude=S)
        if d >= T:
            d = dist(S, N, threads=1) if K == 0 else dist(normalizer(S, N), N, threads=1, exclude=S)
            rec = dict(n=N, k=K, d=d, table=[lo, hi], family="unitary-double-circulant", l=l, g=g, u=u,
                       stab=[list(v) for v in S], time=time.time())
            with open(ROOT / "results" / "qcandidates.jsonl", "a") as f:
                f.write(json.dumps(rec) + "\n")
            print(f"*** QUANTUM CANDIDATE [[{N},{K},{d}]] (table {lo}-{hi}) unitary DC l={l}", flush=True)
    print(f"[[{N},{K}]] unitary-DC: tried {tried}, {time.time() - t0:.0f}s, divisors {len(divs)}",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
