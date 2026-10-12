#!/usr/bin/env python3
"""Which table entries does a new code [n,k,d]_q improve?

Standard propagation rules (all valid for linear codes):
  shortening   [n,k,d] -> [n-1,k-1,d]
  puncturing   [n,k,d] -> [n-1,k,d-1]
  subcode      [n,k,d] -> [n,k-1,d]
  lengthening  [n,k,d] -> [n+1,k,d] (zero column)
Returns the list of (n', k', d', old_lower, upper) with d' > old_lower.
usage: impact.py q n k d
"""
import sys

from qt import load_bounds


def impact(q, n, k, d, B=None):
    B = B or load_bounds(q)
    best = {}
    stack = [(n, k, d)]
    while stack:
        a, b, c = stack.pop()
        if (a, b) not in B or b < 1 or c < 1 or best.get((a, b), 0) >= c:
            continue
        best[(a, b)] = c
        stack += [(a - 1, b - 1, c), (a - 1, b, c - 1), (a, b - 1, c), (a + 1, b, c)]
    out = []
    for (a, b), c in sorted(best.items()):
        lo, hi = B[(a, b)]
        if c > hi:
            raise ValueError(f"[{a},{b},{c}]_{q} would violate the upper bound {hi}: code must be wrong")
        if c > lo:
            out.append((a, b, c, lo, hi))
    return out


if __name__ == "__main__":
    q, n, k, d = map(int, sys.argv[1:5])
    res = impact(q, n, k, d)
    for a, b, c, lo, hi in res:
        print(f"[{a},{b},{c}]_{q}  (table {lo}-{hi})")
    print(f"{len(res)} table entries improved")
