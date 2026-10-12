"""Construction X with an MDS auxiliary code.

Base code C = [n,k,d] (systematic S, pivots piv).  Choose V (k x r) so that no
word of weight < T lies in {u : uV = 0}; append V*A where A generates an
[a, r, delta] code with delta = T - d.  Result: [n+a, k, >= T].
"""
import numpy as np

from gf import GF, avoid, matmul, words_range


def aux_code(q, r, delta):
    """Generator (r x a) of an [a, r, delta]_q code, a = r + delta - 1 (MDS) when possible."""
    F = GF(q)
    a = r + delta - 1
    if delta == 1:
        return np.eye(r, dtype=np.uint8)
    if r == 1:
        return np.ones((1, delta), np.uint8)              # repetition code
    if delta == 2:
        return np.hstack([np.eye(r, dtype=np.uint8), np.ones((r, 1), np.uint8)])  # parity check
    if a > q + 1:
        return None
    # (doubly) extended Reed-Solomon: evaluate x^i at points (+ infinity)
    pts = list(range(q))[: min(a, q)]
    cols = []
    for x in pts:
        col, v = [], 1
        for _ in range(r):
            col.append(v)
            v = int(F.mul[v, x]) if x else 0
        col[0] = 1
        cols.append(col)
    if a == q + 1:
        cols.append([0] * (r - 1) + [1])
    return np.array(cols, np.uint8).T


def options(q, n, k, d, B, tmax=3, rmax=4):
    """Construction-X parameter choices that would beat the table: (T, r, a)."""
    out = []
    for T in range(d + 1, d + tmax + 1):
        for r in range(1, min(rmax, k - 1) + 1):
            A = aux_code(q, r, T - d)
            if A is None:
                continue
            a = A.shape[1]
            lb = B.get((n + a, k))
            if lb and lb[0] < T <= lb[1]:
                out.append((T, r, a))
    return out


def try_x(q, G, d, B, maxwords=200000, threads=2, seed=1, attempts=3000, inner=30000):
    """Yield (T, r, a, Gnew) for every successful Construction X that beats the table."""
    k, n = G.shape
    opts = options(q, n, k, d, B)
    if not opts:
        return
    Tmax = max(T for T, _, _ in opts)
    W, wts, S, piv = words_range(q, G, d, Tmax - 1, maxwords=maxwords, threads=threads)
    if isinstance(W, int):
        # too many words up to Tmax-1: retry with smaller T
        opts = [o for o in opts if o[0] == d + 1]
        if not opts:
            return
        W, wts, S, piv = words_range(q, G, d, d, maxwords=maxwords, threads=threads)
        if isinstance(W, int):
            return
    for T, r, a in sorted(opts, key=lambda o: (o[2], o[1])):
        sel = wts < T
        V = avoid(q, W[sel], piv, k, r, attempts=attempts, inner=inner, seed=seed)
        if V is None:
            continue
        A = aux_code(q, r, T - d)
        Gnew = np.hstack([S, matmul(q, V, A)])
        yield T, r, a, Gnew, int(sel.sum())
