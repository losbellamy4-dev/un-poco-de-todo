#!/usr/bin/env python3
"""Parallel driver for src/qtsearch over many (q, n, k, l, m) targets.

Modes
  survey : short histogram runs (QT_HIST=1) to measure, per target, how often
           random QT codes reach d = lo-2, lo-1, lo (lo = table lower bound).
  hunt   : long fast runs on a ranked target list looking for d >= lo+1.

Every job = (target, lam, divisor g) for `secs` seconds on one core.
Hits are appended to codes/results/candidates.jsonl; per-job stats to --out.

examples
  qt_driver.py survey --q 4 5 8 9 --nmin 60 --maxcost 1e9 --secs 2 --out runs/survey.jsonl
  qt_driver.py hunt --from runs/survey.jsonl --top 40 --secs 120 --out runs/hunt.jsonl
"""
import argparse
import json
import math
import os
import random
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from gf import ROOT
from qt import QTFamily, bz_cost, load_bounds

BIN = str(ROOT / "src" / "qtsearch")
CHAR = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2, 9: 3}
_fam = {}


def family(q, l, lam):
    key = (q, l, lam)
    if key not in _fam:
        _fam[key] = QTFamily(q, l, lam)
    return _fam[key]


def plan_targets(qs, nmin, nmax, kmin, kmax, maxcost, mmax, lmin, flat=False):
    """flat=True: base lengths n with lb(n+1,k) == lb(n,k) < ub(n+1,k) (one-column extension targets)."""
    out = []
    for q in qs:
        B = load_bounds(q)
        for (n, k), (lo, hi) in B.items():
            if not (nmin <= n <= nmax and kmin <= k <= kmax):
                continue
            if flat:
                nxt = B.get((n + 1, k))
                if not nxt or nxt[0] != lo or nxt[1] <= lo:
                    continue
            elif lo >= hi:
                continue
            cost = bz_cost(q, n, k, lo + 1)
            if cost > maxcost:
                continue
            for m in range(2, mmax + 1):
                if n % m:
                    continue
                l = n // m
                if l < k or l < lmin or l % CHAR[q] == 0:
                    continue
                lams = [lam for lam in range(1, q) if (l - k) in family(q, l, lam).reach]
                if lams:
                    out.append(dict(q=q, n=n, k=k, l=l, m=m, lo=lo, hi=hi, cost=cost, lams=lams))
    return out


def ext_flags(t):
    """bit 1: one-column extension of a d=lo code beats lb(n+1,k); bit 2: two columns beat lb(n+2,k)."""
    B = load_bounds_cached(t["q"])
    n, k, lo = t["n"], t["k"], t["lo"]
    f = 0
    a = B.get((n + 1, k))
    if a and a[0] < lo + 1 <= a[1]:
        f |= 1
    b = B.get((n + 2, k))
    if b and b[0] < lo + 2 <= b[1]:
        f |= 2
    return f


_bounds = {}


def load_bounds_cached(q):
    if q not in _bounds:
        _bounds[q] = load_bounds(q)
    return _bounds[q]


def run_job(t, lam, g, secs, seed, hist, ext=False):
    args = [BIN, str(t["q"]), str(t["l"]), str(t["m"]), str(lam), str(t["k"]), str(t["lo"]),
            "0", str(seed), str(secs)] + [str(c) for c in g]
    env = dict(os.environ, OMP_NUM_THREADS="1")
    env.pop("QT_HIST", None)
    env.pop("QT_EXT", None)
    if hist:
        env["QT_HIST"] = "1"
    if ext and ext_flags(t):
        env["QT_EXT"] = str(ext_flags(t))
    try:
        out = subprocess.run(args, capture_output=True, text=True, env=env, timeout=3 * secs + 60).stdout
    except subprocess.TimeoutExpired as e:
        out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
    hits, stats = [], {}
    for line in out.splitlines():
        if line.startswith("EXT"):
            kv = dict(p.split("=") for p in line.split()[1:])
            G = [[int(c) for c in row] for row in kv["matrix"].split(",")]
            hits.append((int(kv["d>"]), G))
        elif line.startswith("HIT"):
            parts = line.split()
            d = int(parts[1][3:])
            fs = [[1]] + [[int(x) for x in p.split("=")[1].split(",")] for p in parts[2:]]
            hits.append((d, fs))
        elif line.startswith("STATS"):
            stats = {kv.split("=")[0]: float(kv.split("=")[1]) for kv in line.split()[1:]}
    return hits, stats


def jobs_for(t, n_jobs, rng):
    jobs = []
    for _ in range(n_jobs):
        lam = rng.choice(t["lams"])
        g = family(t["q"], t["l"], lam).random_divisor(t["l"] - t["k"], rng)
        jobs.append((t, lam, g))
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["survey", "hunt", "ext"])
    ap.add_argument("--q", type=int, nargs="+", default=[7])
    ap.add_argument("--nmin", type=int, default=1)
    ap.add_argument("--nmax", type=int, default=999)
    ap.add_argument("--kmin", type=int, default=3)
    ap.add_argument("--kmax", type=int, default=999)
    ap.add_argument("--lmin", type=int, default=1)
    ap.add_argument("--mmax", type=int, default=8)
    ap.add_argument("--maxcost", type=float, default=1e9)
    ap.add_argument("--secs", type=float, default=2)
    ap.add_argument("--jobs-per-target", type=int, default=1)
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--from", dest="src", default=None, help="survey stats file to rank targets")
    ap.add_argument("--targets", nargs="+", default=None, help="explicit targets q,n,k,l,m")
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=int(time.time()))
    args = ap.parse_args()
    rng = random.Random(args.seed)

    if args.targets:
        targets = []
        for spec in args.targets:
            q, n, k, l, m = map(int, spec.split(","))
            lo, hi = load_bounds_cached(q)[(n, k)]
            lams = [lam for lam in range(1, q) if (l - k) in family(q, l, lam).reach]
            if not lams:
                print(f"skip {spec}: no divisor of degree l-k", file=sys.stderr)
                continue
            targets.append(dict(q=q, n=n, k=k, l=l, m=m, lo=lo, hi=hi, cost=bz_cost(q, n, k, lo + 1), lams=lams))
    elif args.src:
        targets = rank_from_survey(args.src, args.top)
    else:
        targets = plan_targets(args.q, args.nmin, args.nmax, args.kmin, args.kmax,
                               args.maxcost, args.mmax, args.lmin, flat=args.mode == "ext")
    print(f"{len(targets)} targets", file=sys.stderr, flush=True)
    jobs = []
    for t in targets:
        jobs += jobs_for(t, args.jobs_per_target, rng)
    rng.shuffle(jobs)
    cand_path = ROOT / "results" / "candidates.jsonl"
    hist = args.mode == "survey"
    done = 0
    with ThreadPoolExecutor(args.workers) as ex, open(args.out, "a") as fo:
        futs = {ex.submit(run_job, t, lam, g, args.secs, rng.randrange(1 << 30), hist,
                          args.mode == "ext"): (t, lam, g)
                for t, lam, g in jobs}
        for fu in as_completed(futs):
            t, lam, g = futs[fu]
            hits, st = fu.result()
            rec = {k: t[k] for k in ("q", "n", "k", "l", "m", "lo", "hi", "cost")}
            rec.update(lam=lam, g=g, **st)
            fo.write(json.dumps(rec) + "\n")
            fo.flush()
            for d, fs in hits:
                if isinstance(fs[0], list) and len(fs[0]) > t["k"]:   # EXT: full matrix of [n+e,k]
                    c = dict(q=t["q"], n=len(fs[0]), k=t["k"], d_lower_bound=d, kind="matrix", G=fs,
                             origin=dict(kind="QT+ext", n=t["n"], l=t["l"], m=t["m"], lam=lam, g=g),
                             time=time.time())
                    print(f"*** EXT CANDIDATE q={t['q']} [{len(fs[0])},{t['k']}] d>={d}", flush=True)
                    with open(cand_path, "a") as f:
                        f.write(json.dumps(c) + "\n")
                    continue
                c = dict(q=t["q"], n=t["n"], k=t["k"], d_lower_bound=d, table=[t["lo"], t["hi"]],
                         kind="QT", l=t["l"], m=t["m"], lam=lam, g=g, fs=fs, time=time.time())
                with open(cand_path, "a") as f:
                    f.write(json.dumps(c) + "\n")
                print(f"*** CANDIDATE q={t['q']} [{t['n']},{t['k']}] d>={d} table {t['lo']}-{t['hi']}",
                      flush=True)
            done += 1
            if done % 50 == 0:
                print(f"{done}/{len(jobs)} jobs", file=sys.stderr, flush=True)


def rank_from_survey(path, top):
    """Aggregate survey stats per (q,n,k,l,m); score = estimated P(d >= lo+1) per second."""
    agg = {}
    for line in open(path):
        r = json.loads(line)
        if "tried" not in r:
            continue
        key = (r["q"], r["n"], r["k"], r["l"], r["m"])
        a = agg.setdefault(key, dict(r, tried=0, m2=0, m1=0, eq=0, hit=0, secs=0, lams=set()))
        for f in ("tried", "m2", "m1", "eq", "hit", "secs"):
            a[f] += r.get(f, 0)
        a["lams"].add(r["lam"])
    scored = []
    for key, a in agg.items():
        T = max(a["tried"], 1)
        p2, p1, p0 = a["m2"] / T, a["m1"] / T, a["eq"] / T
        if p0 == 0:
            continue
        r1 = p0 / max(p1, 1 / T)
        r0 = p1 / max(p2, 1 / T)
        p_next = p0 * r1 * min(1.0, r1 / max(r0, 1e-9))  # log-concave tail extrapolation
        rate = T / max(a["secs"], 1e-9)
        a["score"] = p_next * rate
        a["p_next"] = p_next
        a["p0"] = p0
        scored.append(a)
    scored.sort(key=lambda a: -a["score"])
    out = []
    for a in scored[:top]:
        q, l, k = a["q"], a["l"], a["k"]
        lams = [lam for lam in range(1, q) if (l - k) in family(q, l, lam).reach]
        out.append(dict(q=q, n=a["n"], k=k, l=l, m=a["m"], lo=a["lo"], hi=a["hi"], cost=a["cost"],
                        lams=lams, score=a["score"]))
        print(f"q={q} [{a['n']},{k}] l={l} m={a['m']} table {a['lo']}-{a['hi']} "
              f"p(lo)={a['p0']:.2e} p_next~{a['p_next']:.1e} rate={a['tried']/max(a['secs'],1e-9):.0f}/s",
              file=sys.stderr)
    return out


if __name__ == "__main__":
    main()
