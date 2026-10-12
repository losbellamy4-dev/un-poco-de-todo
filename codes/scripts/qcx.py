#!/usr/bin/env python3
"""Quantum Construction X (Lisonek-Singh) from nearly Hermitian self-orthogonal GF(4) codes.

C = C0 + X with C0 a Hermitian self-orthogonal QC code [n, kappa-e] (see qqc.py)
and X = e random vectors of C0^perp_h.  If the Hermitian Gram matrix of C has
rank e, there is an [[n+e, n-2*kappa+e, d]] code with
    d >= min( d(C^perp_h), d(C + C^perp_h) + 1 ).
Targets [[N,K]]: kappa = (N-K)/2, n = N - e for e = 1..emax.
Both distances are computed exactly with libmindist (q=4); candidates are
re-verified later.

usage: qcx.py N K [--secs S] [--emax E] [--seed s]
"""
import argparse
import json
import random
import sys
import time

import numpy as np

from gf import GF, ROOT, min_distance
from qqc import CONJ, qc_generator, self_orth_divisors
from qsearch import load_qbounds

F4 = GF(4)
CONJ_A = np.array(CONJ, np.uint8)


def rref(M):
    M = M.copy()
    k, n = M.shape
    r, piv = 0, []
    for col in range(n):
        p = next((i for i in range(r, k) if M[i, col]), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        M[r] = F4.mul[F4.inv[M[r, col]], M[r]]
        for i in range(k):
            if i != r and M[i, col]:
                M[i] = F4.sub[M[i], F4.mul[M[i, col], M[r]]]
        piv.append(col)
        r += 1
        if r == k:
            break
    return M[:r], piv


def nullspace(M, n):
    """Basis of {y : M y^T = 0} (rows)."""
    R, piv = rref(M) if len(M) else (np.zeros((0, n), np.uint8), [])
    free = [c for c in range(n) if c not in piv]
    out = np.zeros((len(free), n), np.uint8)
    for r, f in enumerate(free):
        out[r, f] = 1
        for i, pc in enumerate(piv):
            out[r, pc] = F4.neg[R[i, f]]
    return out


def herm_dual(G, n):
    """C^perp_h = {y : sum x_i conj(y_i) = 0}  <=>  conj(y) in null(G)  <=>  y in conj(null(G))."""
    return CONJ_A[nullspace(G, n)]


def herm_gram_rank(G):
    k, n = G.shape
    Gc = CONJ_A[G]
    gram = np.zeros((k, k), np.uint8)
    for i in range(n):
        gram = F4.add[gram, F4.mul[G[:, i][:, None], Gc[:, i][None, :]]]
    return F4.rank(gram)


def options_for(N, K, emax):
    kappa = (N - K) // 2
    if (N - K) % 2:
        return []
    out = []
    for e in range(1, emax + 1):
        n = N - e
        for m in range(2, 9):
            if n % m:
                continue
            l = n // m
            if l % 2 == 0 or l < 3:
                continue
            try:
                divs = self_orth_divisors(l)
            except ValueError:
                continue
            k0 = kappa - e
            gs = divs.get(l - k0, []) if k0 > 0 else []
            if gs:
                out.append((e, n, m, l, gs))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("N", type=int)
    ap.add_argument("K", type=int)
    ap.add_argument("--secs", type=float, default=300)
    ap.add_argument("--emax", type=int, default=6)
    ap.add_argument("--seed", type=int, default=int(time.time()))
    args = ap.parse_args()
    N, K = args.N, args.K
    B = load_qbounds()
    lo, hi = B[(N, K)]
    T = lo + 1
    kappa = (N - K) // 2
    opts = options_for(N, K, args.emax)
    if not opts:
        print(f"[[{N},{K}]]: no option", file=sys.stderr)
        return
    rng = random.Random(args.seed)
    t0 = time.time()
    tried = best = 0
    out = ROOT / "results" / "qcandidates.jsonl"
    while time.time() - t0 < args.secs:
        e, n, m, l, gs = rng.choice(opts)
        g = rng.choice(gs)
        fs = [[1]] + [[rng.randrange(4) for _ in range(l)] for _ in range(m - 1)]
        G0 = qc_generator(g, fs, l)
        D0 = herm_dual(G0, n)                      # C0^perp_h
        coef = np.array([[rng.randrange(4) for _ in range(len(D0))] for _ in range(e)], np.uint8)
        X = np.zeros((e, n), np.uint8)
        for i in range(len(D0)):
            X = F4.add[X, F4.mul[coef[:, i][:, None], D0[i][None, :]]]
        G = np.vstack([G0, X])
        if F4.rank(G) != kappa or herm_gram_rank(G) != e:
            continue
        tried += 1
        Cd = herm_dual(G, n)                       # [n, n-kappa]
        d1 = min_distance(4, Cd, target=T, threads=1)
        if d1 < T:
            best = max(best, min(d1, T - 1))
            continue
        S = np.vstack([G, Cd])                     # C + C^perp_h
        S, _ = rref(S)
        d2 = min_distance(4, S, target=T - 1, threads=1)
        if d2 + 1 < T:
            continue
        rec = dict(N=N, K=K, n=n, k=K, e=e, d_bound=min(d1, d2 + 1), table=[lo, hi], family="QCX-herm",
                   m=m, l=l, g=g, fs=fs, G=G.tolist(), time=time.time())
        rec["n"] = N
        with open(out, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"*** QUANTUM CANDIDATE [[{N},{K},>={min(d1, d2 + 1)}]] (table {lo}-{hi}) via QCX e={e} "
              f"m={m} l={l}", flush=True)
    print(f"[[{N},{K}]] QCX: tried {tried}, {time.time()-t0:.0f}s, opts {[(e, n, m, l) for e, n, m, l, _ in opts]}",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
