"""Finite fields GF(q), q in {2,3,4,5,7,8,9}, and a ctypes wrapper of libmindist.

Elements are canonical integers a0 + p*a1 + p^2*a2 (coefficients w.r.t. the
Conway polynomial), exactly as in src/mindist.c.
"""
import ctypes
import itertools
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_lib = ctypes.CDLL(str(ROOT / "src" / "libmindist.so"))
_lib.md_compute.restype = ctypes.c_int
_lib.md_compute.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_int, ctypes.c_int,
                            ctypes.c_void_p, ctypes.c_void_p]
_lib.md_field_tables.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]

CHAR = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2, 9: 3}


class GF:
    _cache = {}

    def __new__(cls, q):
        if q in cls._cache:
            return cls._cache[q]
        self = super().__new__(cls)
        self.q, self.p = q, CHAR[q]
        add = np.zeros(q * q, np.uint8)
        mul = np.zeros(q * q, np.uint8)
        inv = np.zeros(q, np.uint8)
        if _lib.md_field_tables(q, add.ctypes.data, mul.ctypes.data, inv.ctypes.data):
            raise ValueError(f"unsupported q={q}")
        self.add, self.mul, self.inv = add.reshape(q, q), mul.reshape(q, q), inv
        self.neg = np.array([int(np.where(self.add[a] == 0)[0][0]) for a in range(q)], np.uint8)
        self.sub = self.add[:, self.neg]  # sub[a,b] = a - b
        # a primitive element and discrete logs (for Magma-compatible output)
        for g in range(2, q) if q > 2 else [1]:
            x, seen = 1, []
            for _ in range(q - 1):
                seen.append(x)
                x = int(self.mul[x, g])
            if len(set(seen)) == q - 1:
                self.prim = g
                self.log = {e: i for i, e in enumerate(seen)}
                self.exp = seen
                break
        cls._cache[q] = self
        return self

    # --- vector helpers (numpy arrays of canonical elements) ---
    def vadd(self, a, b):
        return self.add[a, b]

    def vscale(self, c, a):
        return self.mul[c, a]

    def rank(self, M):
        M = np.array(M, np.uint8)
        k, n = M.shape
        r = 0
        for col in range(n):
            piv = next((i for i in range(r, k) if M[i, col]), None)
            if piv is None:
                continue
            M[[r, piv]] = M[[piv, r]]
            M[r] = self.mul[self.inv[M[r, col]], M[r]]
            for i in range(k):
                if i != r and M[i, col]:
                    M[i] = self.sub[M[i], self.mul[M[i, col], M[r]]]
            r += 1
            if r == k:
                break
        return r


def min_distance(q, G, target=0, threads=0, info=False):
    """Exact minimum distance (target=0), or a test of d >= target (see mindist.c)."""
    G = np.ascontiguousarray(G, dtype=np.uint8)
    k, n = G.shape
    wit = np.zeros(n, np.uint8)
    inf = np.zeros(4, np.int64)
    d = _lib.md_compute(q, n, k, G.ctypes.data, target, threads, wit.ctypes.data, inf.ctypes.data)
    if info:
        return d, wit, dict(enumerated=int(inf[0]), rank=int(inf[1]), infosets=int(inf[2]), w=int(inf[3]))
    return d


def brute_min_distance(q, G):
    """Reference implementation: enumerate all q^k - 1 nonzero codewords."""
    F = GF(q)
    G = np.array(G, np.uint8)
    k, n = G.shape
    best = n + 1
    for coeffs in itertools.product(range(q), repeat=k):
        if not any(coeffs):
            continue
        v = np.zeros(n, np.uint8)
        for c, row in zip(coeffs, G):
            if c:
                v = F.add[v, F.mul[c, row]]
        w = int(np.count_nonzero(v))
        if 0 < w < best:
            best = w
    return best if best <= n else 0


_lib.md_minwords.restype = ctypes.c_int
_lib.md_minwords.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                             ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
_lib.md_extend.restype = ctypes.c_int
_lib.md_extend.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                           ctypes.c_void_p, ctypes.c_double, ctypes.c_longlong, ctypes.c_ulonglong,
                           ctypes.c_void_p]


def min_words(q, G, d, maxwords=100000, threads=0):
    """All weight-d words (projective reps). Returns (words, systematic G, pivots) or (status, None, None)."""
    G = np.ascontiguousarray(G, dtype=np.uint8)
    k, n = G.shape
    W = np.zeros((maxwords, n), np.uint8)
    S = np.zeros((k, n), np.uint8)
    piv = np.zeros(k, np.int32)
    r = _lib.md_minwords(q, n, k, G.ctypes.data, d, maxwords, threads, W.ctypes.data, S.ctypes.data,
                         piv.ctypes.data)
    if r < 0:
        return r, None, None
    return W[:r].copy(), S, piv


def extend_column(q, words, piv, k, exhaustive_limit=5e7, tries=10**7, seed=1):
    """Column v with u.v != 0 for every listed word's message u (or None)."""
    words = np.ascontiguousarray(words, dtype=np.uint8)
    nw, n = words.shape
    piv = np.ascontiguousarray(piv, dtype=np.int32)
    v = np.zeros(k, np.uint8)
    ok = _lib.md_extend(q, n, k, words.ctypes.data, nw, piv.ctypes.data, exhaustive_limit, tries, seed,
                        v.ctypes.data)
    return v if ok == 1 else None


_lib.md_words_range.restype = ctypes.c_int
_lib.md_words_range.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_int,
                                ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p,
                                ctypes.c_void_p, ctypes.c_void_p]
_lib.md_extend2.restype = ctypes.c_int
_lib.md_extend2.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p,
                            ctypes.c_int, ctypes.c_void_p, ctypes.c_int, ctypes.c_double, ctypes.c_longlong,
                            ctypes.c_int, ctypes.c_ulonglong, ctypes.c_void_p, ctypes.c_void_p]


def words_range(q, G, d, dhi, maxwords=200000, threads=0):
    """Projective words with weight in [d, dhi] (d must be the minimum distance)."""
    G = np.ascontiguousarray(G, dtype=np.uint8)
    k, n = G.shape
    W = np.zeros((maxwords, n), np.uint8)
    wts = np.zeros(maxwords, np.int32)
    S = np.zeros((k, n), np.uint8)
    piv = np.zeros(k, np.int32)
    r = _lib.md_words_range(q, n, k, G.ctypes.data, d, dhi, maxwords, threads, W.ctypes.data,
                            wts.ctypes.data, S.ctypes.data, piv.ctypes.data)
    if r < 0:
        return r, None, None, None
    return W[:r].copy(), wts[:r].copy(), S, piv


def extend2(q, words, wts, piv, k, d, exhaustive_limit=2e7, tries=2 * 10**6, poolmax=4000, seed=1):
    words = np.ascontiguousarray(words, dtype=np.uint8)
    wts = np.ascontiguousarray(wts, dtype=np.int32)
    nw, n = words.shape
    piv = np.ascontiguousarray(piv, dtype=np.int32)
    v1 = np.zeros(k, np.uint8)
    v2 = np.zeros(k, np.uint8)
    r = _lib.md_extend2(q, n, k, words.ctypes.data, wts.ctypes.data, nw, piv.ctypes.data, d,
                        exhaustive_limit, tries, poolmax, seed, v1.ctypes.data, v2.ctypes.data)
    return (v1, v2) if r == 1 else None
