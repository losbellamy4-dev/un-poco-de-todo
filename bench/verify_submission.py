#!/usr/bin/env python3
"""Verify a candidate solution for the open-codes benchmark.

Classical submission (JSON):
  {"q": 7, "n": 59, "k": 8, "G": ["1000000065...", ...]}      rows as digit strings
  or "G": [[1,0,0,...], ...].  Field elements are canonical integers
  a0 + p*a1 + p^2*a2  <->  a0 + a1*w + a2*w^2,  w a root of the Conway polynomial
  (x^2+x+1 for GF(4), x^3+x+1 for GF(8), x^2+2x+2 for GF(9)), as in Magma/GAP Z(q).
Quantum submission (JSON):
  {"n": 40, "k": 8, "stabilizer": ["XZZX...", ...]}             n - k Pauli strings over I,X,Y,Z

Output: JSON verdict with the exact minimum distance, the table record and
whether the submission improves it.  --gap additionally recomputes the
distance with GAP/GUAVA (classical codes, or the binary image of a quantum code).
usage: verify_submission.py submission.json [--gap] [--threads N]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "codes" / "scripts"))
from gf import GF, min_distance  # noqa: E402
from qt import load_bounds  # noqa: E402


def load_qbounds():
    B = {}
    for line in open(ROOT.parent / "codes" / "data" / "qbounds_q4.tsv").readlines()[1:]:
        n, k, lo, hi = map(int, line.split())
        B[(n, k)] = (lo, hi)
    return B


def pauli(s):
    x = z = 0
    for i, c in enumerate(s.strip()):
        if c in "XY":
            x |= 1 << i
        if c in "ZY":
            z |= 1 << i
        if c not in "IXYZ":
            raise ValueError(f"bad Pauli character {c!r}")
    return x, z


def verify_classical(sub, use_gap, threads):
    q, n, k = sub["q"], sub["n"], sub["k"]
    rows = sub["G"]
    G = np.array([[int(c) for c in r] if isinstance(r, str) else r for r in rows], np.uint8)
    if G.shape != (k, n):
        return dict(ok=False, error=f"G has shape {G.shape}, expected {(k, n)}")
    if G.max() >= q:
        return dict(ok=False, error="entries must be field elements 0..q-1")
    if GF(q).rank(G) != k:
        return dict(ok=False, error="G does not have full rank k")
    d, wit, info = min_distance(q, G, threads=threads, info=True)
    lo, hi = load_bounds(q)[(n, k)]
    res = dict(ok=True, kind="classical", q=q, n=n, k=k, d=d, record=lo, upper=hi, improves=d > lo,
               witness_weight=int(np.count_nonzero(wit)), enumerated=info["enumerated"])
    if use_gap:
        from verify import gap_min_distance
        gd = gap_min_distance(q, G)
        res["gap"] = dict(dim=gd[0], d=gd[1], agrees=gd[0] == k and gd[1] == d)
    return res


def verify_quantum(sub, use_gap, threads):
    from quantum import normalizer, qdist, rank, sym
    n, k = sub["n"], sub["k"]
    if n > 64:
        return dict(ok=False, error="this verifier supports n <= 64")
    S = [pauli(s) for s in sub["stabilizer"]]
    if any(len(s.strip()) != n for s in sub["stabilizer"]):
        return dict(ok=False, error="every Pauli string must have length n")
    if any(sym(a, b) for a in S for b in S):
        return dict(ok=False, error="stabilizer generators do not commute")
    if rank(S, n) != n - k:
        return dict(ok=False, error=f"stabilizer rank {rank(S, n)} != n - k")
    if k == 0:
        d = qdist(S, n, threads=threads)
        N = S
    else:
        N = normalizer(S, n)
        d = qdist(N, n, threads=threads, exclude=S)
    lo, hi = load_qbounds()[(n, k)]
    res = dict(ok=True, kind="quantum", n=n, k=k, d=d, record=lo, upper=hi, improves=d > lo)
    if use_gap:
        from qverify import binary_image
        from verify import gap_min_distance
        gd = gap_min_distance(2, binary_image(N, n))
        res["gap"] = dict(d_binary_image=gd[1], pure_distance=gd[1] // 2,
                          agrees=(gd[1] == 2 * d) if k == 0 else (gd[1] <= 2 * d))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("submission")
    ap.add_argument("--gap", action="store_true")
    ap.add_argument("--threads", type=int, default=0)
    args = ap.parse_args()
    sub = json.load(open(args.submission))
    try:
        res = verify_quantum(sub, args.gap, args.threads) if "stabilizer" in sub else \
            verify_classical(sub, args.gap, args.threads)
    except Exception as e:  # malformed submissions
        res = dict(ok=False, error=f"{type(e).__name__}: {e}")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
