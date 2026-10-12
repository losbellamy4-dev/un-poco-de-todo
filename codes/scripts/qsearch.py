#!/usr/bin/env python3
"""Search [[n,0,d]] graph codes (circulant / bi-circulant graphs) beating the qubit table.

Each sample is tested against T = lb(n,0)+1 (fast rejection), codes reaching the
record are counted, and every code with d >= T is re-verified exactly and
appended to results/qcandidates.jsonl.

usage: qsearch.py n [--family circ|bicirc] [--secs S] [--seed s] [--target T]
"""
import argparse
import json
import random
import sys
import time

from gf import ROOT
from quantum import circulant_graph, graph_code, qdist


def load_qbounds():
    B = {}
    for line in open(ROOT / "data" / "qbounds_q4.tsv").readlines()[1:]:
        n, k, lo, hi = map(int, line.split())
        B[(n, k)] = (lo, hi)
    return B


def bicirculant_graph(n, A, Bc, C):
    """Two orbits of size h=n/2: adjacency [[circ(A), circ(B)], [circ(B)^T, circ(C)]]."""
    h = n // 2
    adj = [0] * n
    for i in range(h):
        for a in A:
            adj[i] |= 1 << ((i + a) % h)
            adj[i] |= 1 << ((i - a) % h)
        for c in C:
            adj[h + i] |= 1 << (h + (i + c) % h)
            adj[h + i] |= 1 << (h + (i - c) % h)
        for b in Bc:
            j = (i + b) % h
            adj[i] |= 1 << (h + j)
            adj[h + j] |= 1 << i
    for i in range(n):
        adj[i] &= ~(1 << i)
    return adj


def sample(n, family, rng):
    if family == "circ":
        conn = [c for c in range(1, n // 2 + 1) if rng.random() < 0.5]
        return dict(conn=conn), circulant_graph(n, conn)
    h = n // 2
    A = [c for c in range(1, h // 2 + 1) if rng.random() < 0.5]
    C = [c for c in range(1, h // 2 + 1) if rng.random() < 0.5]
    Bc = [c for c in range(h) if rng.random() < 0.5]
    return dict(A=A, B=Bc, C=C), bicirculant_graph(n, A, Bc, C)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("--family", default="circ")
    ap.add_argument("--secs", type=float, default=600)
    ap.add_argument("--seed", type=int, default=int(time.time()))
    ap.add_argument("--target", type=int, default=0)
    args = ap.parse_args()
    n = args.n
    B = load_qbounds()
    lo, hi = B[(n, 0)]
    T = args.target or lo + 1
    rng = random.Random(args.seed)
    t0 = time.time()
    tried = matched = hits = 0
    out = ROOT / "results" / "qcandidates.jsonl"
    while time.time() - t0 < args.secs:
        params, adj = sample(n, args.family, rng)
        g = graph_code(adj, n)
        tried += 1
        d = qdist(g, n, target=T, threads=1)
        if d < T:
            if d >= lo:
                matched += 1
            continue
        d = qdist(g, n, threads=1)
        hits += 1
        rec = dict(n=n, k=0, d=d, table=[lo, hi], family=args.family, params=params, adj=adj, time=time.time())
        with open(out, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"*** QUANTUM CANDIDATE [[{n},0,{d}]] (table {lo}-{hi}) {args.family} {params}", flush=True)
    print(f"n={n} {args.family}: tried {tried}, matched(d>=lo, approx) {matched}, hits {hits}, "
          f"{time.time()-t0:.0f}s", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
