"""Stabilizer (additive GF(4)) codes in symplectic form, n <= 64.

A Pauli is a pair of ints (x, z) (bit i = qubit i).  A stabilizer group is a
list of commuting, independent generators.  Distances via libmindist's
qmd_distance (Brouwer-Zimmermann over qubit positions).
"""
import ctypes
import itertools

import numpy as np

from gf import _lib

_lib.qmd_distance.restype = ctypes.c_int
_lib.qmd_distance.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int,
                              ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_void_p,
                              ctypes.c_void_p]


def sym(a, b):
    """Symplectic inner product of Paulis a=(x,z), b=(x,z)."""
    return (bin(a[0] & b[1]).count("1") + bin(a[1] & b[0]).count("1")) & 1


def weight(v):
    return bin(v[0] | v[1]).count("1")


def rank(gens, n):
    rows = [(x | (z << n)) for x, z in gens]
    r = 0
    for b in range(2 * n):
        p = next((i for i in range(r, len(rows)) if rows[i] >> b & 1), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> b & 1:
                rows[i] ^= rows[r]
        r += 1
    return r


def normalizer(gens, n):
    """Basis of {v : sym(v, g) = 0 for all g} (dimension 2n - rank)."""
    # linear system over GF(2): for each g, sum_i v.x_i g.z_i + v.z_i g.x_i = 0
    m = len(gens)
    rows = []
    for x, z in gens:
        # constraint vector on 2n unknowns (v.x bits then v.z bits): coefficient of v.x_i is g.z_i
        rows.append(z | (x << n))
    # null space of rows (m x 2n)
    piv = []
    R = list(rows)
    r = 0
    for b in range(2 * n):
        p = next((i for i in range(r, m) if R[i] >> b & 1), None)
        if p is None:
            continue
        R[r], R[p] = R[p], R[r]
        for i in range(m):
            if i != r and R[i] >> b & 1:
                R[i] ^= R[r]
        piv.append(b)
        r += 1
    free = [b for b in range(2 * n) if b not in piv]
    basis = []
    for f in free:
        v = 1 << f
        for i, pb in enumerate(piv):
            if R[i] >> f & 1:
                v |= 1 << pb
        basis.append((v & ((1 << n) - 1), v >> n))
    return basis


def qdist(gens, n, target=0, threads=0, exclude=None, info=False):
    """Min weight of span(gens) \\ span(exclude)."""
    X = np.array([g[0] for g in gens], np.uint64)
    Z = np.array([g[1] for g in gens], np.uint64)
    ex = exclude or []
    SX = np.array([g[0] for g in ex] or [0], np.uint64)
    SZ = np.array([g[1] for g in ex] or [0], np.uint64)
    wit = np.zeros(2, np.uint64)
    inf = np.zeros(4, np.int64)
    d = _lib.qmd_distance(n, len(gens), X.ctypes.data, Z.ctypes.data, len(ex), SX.ctypes.data, SZ.ctypes.data,
                          target, threads, wit.ctypes.data, inf.ctypes.data)
    if info:
        return d, (int(wit[0]), int(wit[1])), dict(enumerated=int(inf[0]), rank=int(inf[1]), sets=int(inf[2]),
                                                    w=int(inf[3]))
    return d


def quantum_distance(stab, n, target=0, threads=0):
    """Distance of the stabilizer code [[n, n - rank(stab)]]: min wt(N(S) \\ S) (or min wt(S) if k = 0)."""
    k = n - rank(stab, n)
    if k == 0:
        return qdist(stab, n, target, threads)
    return qdist(normalizer(stab, n), n, target, threads, exclude=stab)


def brute_quantum_distance(stab, n):
    k = n - rank(stab, n)
    N = stab if k == 0 else normalizer(stab, n)
    S = set()
    for coeffs in itertools.product([0, 1], repeat=len(stab)):
        x = z = 0
        for c, g in zip(coeffs, stab):
            if c:
                x ^= g[0]
                z ^= g[1]
        S.add((x, z))
    best = n + 1
    for coeffs in itertools.product([0, 1], repeat=len(N)):
        x = z = 0
        for c, g in zip(coeffs, N):
            if c:
                x ^= g[0]
                z ^= g[1]
        if (x, z) == (0, 0) or (k > 0 and (x, z) in S):
            continue
        best = min(best, bin(x | z).count("1"))
    return best if best <= n else 0


def graph_code(adj, n):
    """[[n,0]] graph-state stabilizer: X_i Z^{N(i)}; adj = list of neighbour bitmasks."""
    return [(1 << i, adj[i]) for i in range(n)]


def circulant_graph(n, conn):
    """Circulant graph on Z_n with connection set conn (symmetrised)."""
    adj = []
    for i in range(n):
        m = 0
        for c in conn:
            m |= 1 << ((i + c) % n)
            m |= 1 << ((i - c) % n)
        m &= ~(1 << i)
        adj.append(m)
    return adj
