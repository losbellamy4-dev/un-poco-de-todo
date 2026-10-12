#!/usr/bin/env python3
"""Try to extend the table's own record codes at flat spots.

For each (n,k) with lb(n+1,k) == lb(n,k) < ub(n+1,k): rebuild the record code C
= [n,k,lo] from codetables.de, list its minimum-weight words and search a
column v with u.v != 0 for all of them  ->  [n+1, k, lo+1] beats the table.
If lb(n+2,k) <= lo+1 also try two columns (-> [n+2,k,lo+2]).
Candidates are appended to results/candidates.jsonl.

usage: extend_table.py q [--maxcost 1e9] [--maxwords 20000] [--nmin N --nmax N]
"""
import argparse
import json
import sys
import time

import numpy as np

from constructions import flat_spots
from gf import ROOT, extend2, extend_column, min_distance, min_words, words_range
from qt import bz_cost, load_bounds
from table_codes import UnsupportedStep, build_entry


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int)
    ap.add_argument("--maxcost", type=float, default=1e9)
    ap.add_argument("--maxwords", type=int, default=20000)
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=999)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--log", default=None)
    args = ap.parse_args()
    q = args.q
    B = load_bounds(q)
    log = open(args.log, "a") if args.log else None
    spots = [(n, k) for n, k in flat_spots(q) if args.nmin <= n <= args.nmax]
    spots = [(n, k) for n, k in spots if bz_cost(q, n, k, B[(n, k)][0] + 1) <= args.maxcost]
    print(f"q={q}: {len(spots)} flat spots within cost", file=sys.stderr, flush=True)
    stats = dict(built=0, unsupported=0, tried1=0, ok1=0, tried2=0, ok2=0)
    cand_path = ROOT / "results" / "candidates.jsonl"
    for i, (n, k) in enumerate(spots):
        lo = B[(n, k)][0]
        rec = dict(q=q, n=n, k=k, lo=lo)
        try:
            G, dclaim = build_entry(q, n, k)
        except (UnsupportedStep, ValueError, KeyError, IndexError) as e:
            stats["unsupported"] += 1
            rec["unsupported"] = str(e)[:100]
            if log:
                log.write(json.dumps(rec) + "\n")
            continue
        stats["built"] += 1
        d = min_distance(q, G, threads=args.threads)
        rec["d"] = d
        if d != lo:
            rec["note"] = f"rebuilt code has d={d}, table says {lo}"
        if d < lo:
            if log:
                log.write(json.dumps(rec) + "\n")
            continue
        W, S, piv = min_words(q, G, d, maxwords=args.maxwords, threads=args.threads)
        rec["A_d"] = len(W) if not isinstance(W, int) else W
        if not isinstance(W, int):
            stats["tried1"] += 1
            v = extend_column(q, W, piv, k, exhaustive_limit=3e7, tries=5 * 10**6, seed=i + 1)
            if v is not None:
                G1 = np.hstack([S, v.reshape(k, 1)])
                stats["ok1"] += 1
                rec["ext1"] = True
                c = dict(q=q, n=n + 1, k=k, d_lower_bound=d + 1, kind="matrix", G=G1.tolist(),
                         origin=dict(kind="table-code+1col", n=n), time=time.time())
                with open(cand_path, "a") as f:
                    f.write(json.dumps(c) + "\n")
                print(f"*** EXT CANDIDATE q={q} [{n+1},{k}] d>={d+1} (table {B.get((n+1,k))})", flush=True)
        nb2 = B.get((n + 2, k))
        if nb2 and nb2[0] < d + 2 <= nb2[1]:
            W2, wts, S2, piv2 = words_range(q, G, d, d + 1, maxwords=10 * args.maxwords, threads=args.threads)
            if not isinstance(W2, int):
                stats["tried2"] += 1
                rec["A_d1"] = int((wts == d + 1).sum())
                r = extend2(q, W2, wts, piv2, k, d, seed=i + 7)
                if r is not None:
                    G2 = np.hstack([S2, r[0].reshape(k, 1), r[1].reshape(k, 1)])
                    stats["ok2"] += 1
                    rec["ext2"] = True
                    c = dict(q=q, n=n + 2, k=k, d_lower_bound=d + 2, kind="matrix", G=G2.tolist(),
                             origin=dict(kind="table-code+2col", n=n), time=time.time())
                    with open(cand_path, "a") as f:
                        f.write(json.dumps(c) + "\n")
                    print(f"*** EXT2 CANDIDATE q={q} [{n+2},{k}] d>={d+2} (table {nb2})", flush=True)
        if log:
            log.write(json.dumps(rec) + "\n")
            log.flush()
        if i % 20 == 0:
            print(f"{i}/{len(spots)} {stats}", file=sys.stderr, flush=True)
    print(f"q={q} done: {stats}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
