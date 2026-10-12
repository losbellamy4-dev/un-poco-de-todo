#!/usr/bin/env python3
"""Cyclic additive stabilizer codes: S = {(c, c*s) : c in <g>} over GF(2)[x]/(x^n - 1).

* s palindromic (s(x) = s(x^-1)) and s_0 = 0  => all elements commute (circulant graph state).
* g | x^n - 1 of degree k  => S has rank n - k, i.e. an [[n, k]] code.
k = 0 (g = 1) gives circulant graph codes.  Even n is handled through repeated
factors: x^n - 1 = (x^m - 1)^(2^e) with m odd.

usage: qcyclic.py n k [--secs S] [--seed s] [--target T]
"""
import argparse
import json
import random
import sys
import time

from gf import ROOT
from poly import PolyRing
from quantum import qdist, normalizer, rank
from qsearch import load_qbounds


def rot(c, j, n):
    j %= n
    mask = (1 << n) - 1
    return ((c << j) | (c >> (n - j))) & mask


def mulmod(a, b, n):
    """Product of binary polynomials (bitmasks) modulo x^n - 1."""
    r = 0
    j = 0
    while b:
        if b & 1:
            r ^= rot(a, j, n)
        b >>= 1
        j += 1
    return r


def factor_degrees(n):
    """Irreducible factors of x^m - 1 over GF(2) (m odd part of n) and multiplicity bound 2^e."""
    e, m = 0, n
    while m % 2 == 0:
        m //= 2
        e += 1
    R = PolyRing(2)
    facs = R.factor_squarefree(R.x_n_minus(m, 1))
    return [sum(c << i for i, c in enumerate(f)) for f in facs], [len(f) - 1 for f in facs], 2 ** e


def pmul2(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        b >>= 1
        a <<= 1
    return r


class CyclicFamily:
    def __init__(self, n):
        self.n = n
        self.facs, self.degs, self.mult = factor_degrees(n)
        reach = {0}
        for d in self.degs:
            reach = {r + d * t for r in reach for t in range(self.mult + 1)}
        self.reach = reach

    def random_g(self, k, rng, tries=2000):
        idx = list(range(len(self.facs)))
        for _ in range(tries):
            rng.shuffle(idx)
            tot, g, chosen = 0, 1, []
            for i in idx:
                for _t in range(self.mult):
                    if tot + self.degs[i] <= k and rng.random() < 0.7:
                        g = pmul2(g, self.facs[i])
                        tot += self.degs[i]
                        chosen.append(i)
                if tot == k:
                    return g, chosen
        return None, None

    def stabilizer(self, g, s):
        n = self.n
        k = g.bit_length() - 1
        gens = []
        for i in range(n - k):
            c = rot(g, i, n)
            gens.append((c, mulmod(c, s, n)))
        return gens


def random_s(n, rng, p=0.5):
    s = 0
    for j in range(1, n // 2 + 1):
        if rng.random() < p:
            s |= 1 << j
            s |= 1 << ((n - j) % n)
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("k", type=int)
    ap.add_argument("--secs", type=float, default=600)
    ap.add_argument("--seed", type=int, default=int(time.time()))
    ap.add_argument("--target", type=int, default=0)
    args = ap.parse_args()
    n, k = args.n, args.k
    B = load_qbounds()
    lo, hi = B[(n, k)]
    T = args.target or lo + 1
    fam = CyclicFamily(n)
    if k not in fam.reach:
        print(f"k={k} not reachable for n={n} (degrees {fam.degs}, mult {fam.mult})", file=sys.stderr)
        return
    rng = random.Random(args.seed)
    t0 = time.time()
    tried = matched = hits = 0
    out = ROOT / "results" / "qcandidates.jsonl"
    while time.time() - t0 < args.secs:
        g, _ = fam.random_g(k, rng)
        if g is None:
            continue
        s = random_s(n, rng)
        S = fam.stabilizer(g, s)
        tried += 1
        if k == 0:
            d = qdist(S, n, target=T, threads=1)
        else:
            N = normalizer(S, n)
            d = qdist(N, n, target=T, threads=1, exclude=S)
        if d < T:
            matched += d >= lo
            continue
        assert rank(S, n) == n - k
        d = qdist(S, n, threads=1) if k == 0 else qdist(normalizer(S, n), n, threads=1, exclude=S)
        hits += 1
        rec = dict(n=n, k=k, d=d, table=[lo, hi], family="cyclic-additive", g=g, s=s, time=time.time())
        with open(out, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"*** QUANTUM CANDIDATE [[{n},{k},{d}]] (table {lo}-{hi}) g={g:#x} s={s:#x}", flush=True)
    print(f"[[{n},{k}]] cyclic: tried {tried}, matched {matched}, hits {hits}, {time.time()-t0:.0f}s",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
