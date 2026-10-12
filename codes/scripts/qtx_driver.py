#!/usr/bin/env python3
"""Driver for src/qtx: local search for QT codes [n=m*l, k] whose j-column
extension (j = 0, 1, 2) beats the table at [n+j, k].

Target selection: base (n, k) and j with T = lb(n+j,k) + 1 <= ub(n+j,k).
slack s = (T - j) - lb(n,k) is how much the base must exceed the table at n:
  s = 1  : base must itself be a record (only j = 0 is kept for these)
  s = 0  : base must match the record at n
  s = -1 : base one below the record at n suffices (j = 2, "double flat")
Targets are ranked by slack, then by verification cost.

usage: qtx_driver.py --q 7 --nmax 70 --maxcost 1e8 --secs 60 --out runs/qtx.jsonl [--slack -1 0]
"""
import argparse
import json
import os
import random
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from gf import ROOT
from qt import QTFamily, bz_cost, load_bounds
from gf import matmul
from xext import aux_code
import numpy as np

BIN = str(ROOT / "src" / "qtX")
CHAR = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2, 9: 3}
_fam = {}


def family(q, l, lam):
    key = (q, l, lam)
    if key not in _fam:
        _fam[key] = QTFamily(q, l, lam)
    return _fam[key]


def plan(qs, nmin, nmax, kmin, kmax, maxcost, mmax, slacks, js, rs=(0,)):
    """js = deltas; rs = subcode codimensions (0 = plain column extension, needs delta <= 2)."""
    out = []
    for q in qs:
        B = load_bounds(q)
        for (n, k), (lo, hi) in B.items():
            if not (nmin <= n <= nmax and kmin <= k <= kmax):
                continue
            for j in js:
                for r in rs:
                    if r == 0:
                        if j > 2:
                            continue
                        a = j
                    else:
                        A = aux_code(q, r, j)
                        if A is None or r >= k:
                            continue
                        a = A.shape[1]
                        if r == 1 and j > 1:
                            continue      # same as plain extension with repetition: weaker
                    nb = B.get((n + a, k))
                    if not nb:
                        continue
                    T = nb[0] + 1
                    if T > nb[1]:
                        continue
                    s_ = (T - j) - lo
                    if s_ not in slacks:
                        continue
                    cost = bz_cost(q, n, k, T)
                    if cost > maxcost:
                        continue
                    for m in range(2, mmax + 1):
                        if n % m:
                            continue
                        l = n // m
                        if l < k or l % CHAR[q] == 0:
                            continue
                        lams = [lam for lam in range(1, q) if (l - k) in family(q, l, lam).reach]
                        if lams:
                            out.append(dict(q=q, n=n, k=k, j=j, r=r, a=a, T=T, s=s_, l=l, m=m, lo=lo, cost=cost,
                                            lams=lams, table_target=list(nb)))
    out.sort(key=lambda t: (t["s"], t["cost"]))
    return out


def run(t, lam, g, secs, seed, patience, ftry):
    args = [BIN, str(t["q"]), str(t["l"]), str(t["m"]), str(lam), str(t["k"]), str(t["T"]), str(t["j"]),
            str(t["r"]), str(secs), str(seed), str(patience), str(ftry)] + [str(c) for c in g]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    try:
        out = subprocess.run(args, capture_output=True, text=True, env=env, timeout=3 * secs + 120).stdout
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--q", type=int, nargs="+", default=[7])
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=999)
    ap.add_argument("--kmin", type=int, default=4)
    ap.add_argument("--kmax", type=int, default=999)
    ap.add_argument("--mmax", type=int, default=8)
    ap.add_argument("--maxcost", type=float, default=1e8)
    ap.add_argument("--slack", type=int, nargs="+", default=[-1, 0])
    ap.add_argument("--j", type=int, nargs="+", default=[1, 2], help="delta = T - base d")
    ap.add_argument("--r", type=int, nargs="+", default=[0], help="subcode codimension (0 = plain extension)")
    ap.add_argument("--secs", type=float, default=60)
    ap.add_argument("--patience", type=int, default=300)
    ap.add_argument("--ftry", type=float, default=-4.0)
    ap.add_argument("--jobs-per-target", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0, help="only the first N targets")
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=int(time.time()))
    ap.add_argument("--dry", action="store_true")
    args = ap.parse_args()
    rng = random.Random(args.seed)
    T = plan(args.q, args.nmin, args.nmax, args.kmin, args.kmax, args.maxcost, args.mmax, set(args.slack),
             args.j, args.r)
    if args.limit:
        T = T[: args.limit]
    print(f"{len(T)} targets", file=sys.stderr, flush=True)
    if args.dry:
        for t in T[:50]:
            print(t, file=sys.stderr)
        return
    jobs = []
    for t in T:
        for _ in range(args.jobs_per_target):
            lam = rng.choice(t["lams"])
            g = family(t["q"], t["l"], lam).random_divisor(t["l"] - t["k"], rng)
            jobs.append((t, lam, g))
    cand = ROOT / "results" / "candidates.jsonl"
    done = 0
    with ThreadPoolExecutor(args.workers) as ex, open(args.out, "a") as fo:
        futs = {ex.submit(run, t, lam, g, args.secs, rng.randrange(1 << 30), args.patience, args.ftry): (t, lam, g)
                for t, lam, g in jobs}
        for fu in as_completed(futs):
            t, lam, g = futs[fu]
            out = fu.result()
            rec = {k: t[k] for k in ("q", "n", "k", "j", "r", "a", "T", "s", "l", "m", "lo", "cost", "table_target")}
            rec.update(lam=lam, g=g)
            best = []
            for line in out.splitlines():
                if line.startswith("STATS"):
                    rec.update({kv.split("=")[0]: float(kv.split("=")[1]) for kv in line.split()[1:]})
                elif line.startswith("BEST"):
                    best.append(line[5:])
                elif line.startswith("XV"):
                    kv = dict(p.split("=", 1) for p in line.split()[1:])
                    S = np.array([[int(c) for c in row] for row in kv["S"].split(",")], np.uint8)
                    V = np.array([[int(c) for c in col] for col in kv["V"].split(",")], np.uint8).T
                    A = aux_code(t["q"], t["r"], t["j"])
                    G = np.hstack([S, matmul(t["q"], V, A)]).tolist()
                    c = dict(q=t["q"], n=len(G[0]), k=t["k"], d_lower_bound=t["T"], kind="matrix", G=G,
                             origin=dict(kind=f"QT+X(r={t['r']},delta={t['j']})", n=t["n"], l=t["l"], m=t["m"],
                                         lam=lam, g=g), time=time.time())
                    with open(cand, "a") as f:
                        f.write(json.dumps(c) + "\n")
                    print(f"*** X CANDIDATE q={t['q']} [{len(G[0])},{t['k']}] d>={t['T']} "
                          f"(table {t['table_target']})", flush=True)
                elif line.startswith("EXT"):
                    kv = dict(p.split("=", 1) for p in line.split()[1:])
                    G = [[int(c) for c in row] for row in kv["matrix"].split(",")]
                    c = dict(q=t["q"], n=len(G[0]), k=t["k"], d_lower_bound=int(kv["d>"].lstrip("=")),
                             kind="matrix", G=G, origin=dict(kind=f"QT+{t['j']}col", n=t["n"], l=t["l"],
                                                             m=t["m"], lam=lam, g=g), time=time.time())
                    with open(cand, "a") as f:
                        f.write(json.dumps(c) + "\n")
                    print(f"*** EXT CANDIDATE q={t['q']} [{len(G[0])},{t['k']}] d>={c['d_lower_bound']} "
                          f"(table {t['table_target']})", flush=True)
                elif line.startswith("HIT"):
                    fs = [[1]] + [[int(x) for x in part.split(",")] for part in line.split("f=")[1].split(";")]
                    c = dict(q=t["q"], n=t["n"], k=t["k"], d_lower_bound=t["T"], kind="QT", l=t["l"], m=t["m"],
                             lam=lam, g=g, fs=fs, time=time.time())
                    with open(cand, "a") as f:
                        f.write(json.dumps(c) + "\n")
                    print(f"*** HIT CANDIDATE q={t['q']} [{t['n']},{t['k']}] d>={t['T']}", flush=True)
            rec["best"] = best[-5:]
            fo.write(json.dumps(rec) + "\n")
            fo.flush()
            done += 1
            if done % 20 == 0:
                print(f"{done}/{len(jobs)} jobs", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
