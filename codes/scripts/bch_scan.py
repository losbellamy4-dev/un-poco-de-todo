#!/usr/bin/env python3
"""Certified-distance scan of constacyclic BCH-type codes against the tables.

lambda in GF(q)* of order r (r | q-1), gcd(n, q) = 1.  The roots of x^n - lambda
are beta^(1 + r*j), j in Z_n, beta a primitive (r*n)-th root of unity; the
q-cyclotomic classes act on J = {j} via j -> (q*(1+r*j) - 1)/r mod n.
A defining set containing {j0 + i*b : 0 <= i < delta-1} with gcd(b, n) = 1
certifies d >= delta (BCH bound).  For every run we take the smallest union of
classes containing it, giving a code [n, n - |Z|, >= delta].
Prints every (n, k) where delta exceeds the table's lower bound.
usage: bch_scan.py q [--nmax N]
"""
import math
import sys
from collections import defaultdict

from qt import load_bounds


def order(a, m):
    if math.gcd(a, m) != 1:
        return None
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k


def classes(n, q, r):
    """Partition of Z_n (indices j of roots beta^(1+rj)) into q-cyclotomic classes."""
    N = r * n
    cls = [-1] * n
    reps = []
    for j in range(n):
        if cls[j] >= 0:
            continue
        cid = len(reps)
        reps.append([])
        e = (1 + r * j) % N
        while True:
            jj = ((e - 1) // r) % n
            if cls[jj] >= 0:
                break
            cls[jj] = cid
            reps[cid].append(jj)
            e = e * q % N
    return cls, [len(c) for c in reps]


def scan(q, nmax):
    B = load_bounds(q)
    p = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2, 9: 3}[q]
    best = defaultdict(int)          # (n, k) -> certified delta
    for n in range(2, nmax + 1):
        if n % p == 0:
            continue
        for r in [d for d in range(1, q) if (q - 1) % d == 0]:
            cls, sizes = classes(n, q, r)
            units = [b for b in range(1, n) if math.gcd(b, n) == 1] or [1]
            seen_b = set()
            for b in units:
                if (n - b) % n in seen_b:      # step -b gives the same runs reversed
                    continue
                seen_b.add(b)
                for j0 in range(n):
                    used = set()
                    size = 0
                    for L in range(1, n):     # run of L consecutive elements -> delta = L + 1
                        j = (j0 + (L - 1) * b) % n
                        c = cls[j]
                        if c not in used:
                            used.add(c)
                            size += sizes[c]
                        k = n - size
                        if k <= 0:
                            break
                        if L + 1 > best[(n, k)]:
                            best[(n, k)] = L + 1
    # propagate: [n,k,d] -> [n,k-1,d]
    improvements = []
    for (n, k), d in sorted(best.items()):
        for kk in range(k, 0, -1):
            if (n, kk) in B and d > B[(n, kk)][0]:
                improvements.append((n, kk, d, B[(n, kk)]))
    return best, improvements


if __name__ == "__main__":
    q = int(sys.argv[1])
    nmax = int(sys.argv[sys.argv.index("--nmax") + 1]) if "--nmax" in sys.argv else 999
    nmax = min(nmax, max(n for n, _ in load_bounds(q)))
    best, imp = scan(q, nmax)
    print(f"q={q}: {len(best)} (n,k) pairs with BCH-certified codes up to n={nmax}; "
          f"{len(imp)} beat the table")
    for t in imp[:30]:
        print("  [%d,%d] certified d>=%d, table %s" % t)
