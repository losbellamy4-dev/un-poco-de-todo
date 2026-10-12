#!/usr/bin/env python3
"""Construction X applied to codetables.de's own record codes.

For every entry (n,k) whose record code C = [n,k,lo] admits a Construction X
option (T, r, a) with T > lb(n+a,k) (see xext.options), rebuild C from its
published construction and search V.  Successes go to results/candidates.jsonl.

usage: xtable.py q [--maxcost 1e9] [--nmin --nmax] [--log runs/x.jsonl] [--shard i/N]
"""
import argparse
import json
import sys
import time

from gf import ROOT, min_distance
from qt import bz_cost, load_bounds
from table_codes import UnsupportedStep, build_entry
from xext import options, try_x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("q", type=int)
    ap.add_argument("--maxcost", type=float, default=1e9)
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=999)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--log", default=None)
    args = ap.parse_args()
    q = args.q
    B = load_bounds(q)
    si, sn = map(int, args.shard.split("/"))
    entries = []
    for (n, k), (lo, hi) in sorted(B.items()):
        if k < 2 or not (args.nmin <= n <= args.nmax):
            continue
        if options(q, n, k, lo, B) and bz_cost(q, n, k, lo + 3) <= args.maxcost:
            entries.append((n, k))
    entries = entries[si::sn]
    print(f"q={q}: {len(entries)} entries", file=sys.stderr, flush=True)
    log = open(args.log, "a") if args.log else None
    st = dict(built=0, unsupported=0, wrong_d=0, success=0)
    cand_path = ROOT / "results" / "candidates.jsonl"
    for i, (n, k) in enumerate(entries):
        lo = B[(n, k)][0]
        rec = dict(q=q, n=n, k=k, lo=lo, opts=options(q, n, k, lo, B))
        try:
            G, _ = build_entry(q, n, k)
        except Exception as e:  # unsupported or unparsable construction
            st["unsupported"] += 1
            rec["unsupported"] = f"{type(e).__name__}: {str(e)[:80]}"
            if log:
                log.write(json.dumps(rec) + "\n")
            continue
        st["built"] += 1
        d = min_distance(q, G, threads=args.threads)
        rec["d"] = d
        if d != lo:
            st["wrong_d"] += 1
        if d >= lo:
            for T, r, a, Gnew, nw in try_x(q, G, d, B, threads=args.threads, seed=i + 1):
                st["success"] += 1
                rec.setdefault("success", []).append([T, r, a, nw])
                c = dict(q=q, n=n + a, k=k, d_lower_bound=T, kind="matrix", G=Gnew.tolist(),
                         origin=dict(kind=f"table-code+X(r={r},a={a})", n=n), time=time.time())
                with open(cand_path, "a") as f:
                    f.write(json.dumps(c) + "\n")
                print(f"*** X CANDIDATE q={q} [{n+a},{k}] d>={T} (table {B[(n+a,k)]}) from [{n},{k},{d}] r={r}",
                      flush=True)
        if log:
            log.write(json.dumps(rec) + "\n")
            log.flush()
        if i % 25 == 0:
            print(f"{i}/{len(entries)} {st}", file=sys.stderr, flush=True)
    print(f"q={q} done: {st}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
