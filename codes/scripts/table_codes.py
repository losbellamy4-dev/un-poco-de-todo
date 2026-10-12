#!/usr/bin/env python3
"""Rebuild codetables.de record codes from their (Magma) construction text.

Supported steps: stored generator matrix, QuasiCyclicCode, QuasiTwistedCyclicCode,
CyclicCode, puncturing, shortening.  Anything else -> UnsupportedStep.
Field elements: integers (prime subfield), w (primitive element = root of the
Conway polynomial, i.e. canonical element p) and products/powers thereof.
"""
import re

import numpy as np

from constructions import fetch, parse
from gf import GF


class UnsupportedStep(Exception):
    pass


def elem(F, tok):
    """Parse a field element like '3', 'w', 'w^5', '2*w^3', '-1'."""
    tok = tok.strip()
    val = 1
    neg = tok.startswith("-")
    if neg:
        tok = tok[1:]
    for f in tok.split("*"):
        f = f.strip()
        if not f:
            continue
        if f.startswith("w"):
            e = int(f[2:]) if f.startswith("w^") else 1
            x = 1
            for _ in range(e):
                x = int(F.mul[x, F.p])        # canonical index of w is p
            val = int(F.mul[val, x])
        elif f.startswith("Z("):
            raise UnsupportedStep("Z() notation")
        else:
            c = int(f) % F.p
            val = int(F.mul[val, c])
    return int(F.neg[val]) if neg else val


def poly(F, s, length):
    """Parse a polynomial string in x into a coefficient vector of given length."""
    coeffs = np.zeros(length, np.uint8)
    s = s.replace(" ", "")
    if not s:
        return coeffs
    # split into signed terms
    terms = re.findall(r"[+-]?[^+-]+", s)
    for t in terms:
        sign = -1 if t.startswith("-") else 1
        t = t.lstrip("+-")
        if "x" in t:
            pre, _, post = t.partition("x")
            e = int(post[1:]) if post.startswith("^") else 1
            c = elem(F, pre.rstrip("*")) if pre.rstrip("*") else 1
        else:
            e, c = 0, elem(F, t)
        if sign < 0:
            c = int(F.neg[c])
        if e >= length:
            raise UnsupportedStep("degree too large")
        coeffs[e] = F.add[coeffs[e], c]
    return coeffs


def twist_rows(F, p, l, lam, rows=None):
    rows = rows or l
    out = np.zeros((rows, l), np.uint8)
    v = p.copy()
    for i in range(rows):
        out[i] = v
        last = v[-1]
        v = np.roll(v, 1)
        v[0] = F.mul[lam, last]
    return out


def row_reduce(F, M):
    """Row echelon form, rank, pivots."""
    M = M.copy()
    k, n = M.shape
    r, piv = 0, []
    for col in range(n):
        p = next((i for i in range(r, k) if M[i, col]), None)
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        M[r] = F.mul[F.inv[M[r, col]], M[r]]
        for i in range(k):
            if i != r and M[i, col]:
                M[i] = F.sub[M[i], F.mul[M[i, col], M[r]]]
        piv.append(col)
        r += 1
        if r == k:
            break
    return M[:r], piv


def positions(s):
    """'{ 111 .. 114 }' or '{ 1, 4, 7 .. 9 }' -> 0-based sorted list."""
    inner = s[s.index("{") + 1: s.index("}")]
    out = []
    for part in inner.split(","):
        part = part.strip()
        if ".." in part:
            a, b = part.split("..")
            out += list(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return sorted(p - 1 for p in out)


def shorten(F, G, pos):
    """Subcode vanishing on pos, then delete pos."""
    k, n = G.shape
    order = list(pos) + [c for c in range(n) if c not in pos]
    M, piv = row_reduce(F, G[:, order])
    # rows with pivot inside the first len(pos) columns are removed
    keep = [i for i, c in enumerate(piv) if c >= len(pos)]
    M = M[keep][:, len(pos):]
    inv_order = [c for c in range(n) if c not in pos]
    return M  # columns already in original order of remaining coordinates


def dual(F, G):
    """Generator matrix of the dual code."""
    M, piv = row_reduce(F, G)
    k, n = M.shape
    free = [c for c in range(n) if c not in piv]
    H = np.zeros((len(free), n), np.uint8)
    for r, f in enumerate(free):
        H[r, f] = 1
        for i, pc in enumerate(piv):
            H[r, pc] = F.neg[M[i, f]]
    return H


def build_step(F, step, built):
    desc, body = step["desc"], " ".join(step["body"])
    n = step["n"]
    if "stored generator matrix" in body:
        rows = re.findall(r"\[([^\[\]]+)\]", body)
        G = np.array([[elem(F, x) for x in r.split(",")] for r in rows], np.uint8)
        return G
    m = re.search(r"QuasiTwistedCyclicCode of length (\d+) and constant (\S+) with generators: (.*)", body)
    if m:
        N, lam = int(m[1]), elem(F, m[2])
        polys = [p for p in m[3].split(",") if p.strip()]
        l = N // len(polys)
        return np.hstack([twist_rows(F, poly(F, p, l), l, lam) for p in polys])
    m = re.search(r"QuasiCyclicCode of length (\d+) with generating polynomials: (.*)", body)
    if m:
        N = int(m[1])
        polys = [p for p in m[2].split(",") if p.strip()]
        l = N // len(polys)
        return np.hstack([twist_rows(F, poly(F, p, l), l, 1) for p in polys])
    m = re.search(r"CyclicCode of length (\d+) with generating polynomial (.*)", body)
    if m:
        N = int(m[1])
        g = poly(F, m[2], N)
        deg = int(np.nonzero(g)[0].max())
        return twist_rows(F, g, N, 1, rows=N - deg)
    m = re.search(r"QuasiCyclicCode of length (\d+) stacked to height (\d+) with generating polynomials: (.*)", body)
    if m:
        N, h = int(m[1]), int(m[2])
        polys = [p for p in m[3].split(",") if p.strip()]
        per = len(polys) // h
        l = N // per
        rows = [np.hstack([twist_rows(F, poly(F, polys[i * per + j], l), l, 1) for j in range(per)])
                for i in range(h)]
        return np.vstack(rows)
    m = re.search(r"ConstaCyclicCode generated by (.*) with shift constant (\S+)", body)
    if m:
        lam = elem(F, m[2])
        g = poly(F, m[1], n)
        deg = int(np.nonzero(g)[0].max())
        return twist_rows(F, g, n, lam, rows=n - deg)
    m = re.search(r"RepetitionCode of length (\d+)", body)
    if m and "Dual" not in body:
        return np.ones((1, int(m[1])), np.uint8)
    m = re.search(r"^Subcode of \[(\d+)\]", body)
    if m:
        G, _ = row_reduce(F, built[int(m[1])])
        return G[: step["k"]]
    m = re.search(r"PadCode \[(\d+)\] by (\d+)", body)
    if m:
        G = built[int(m[1])]
        return np.hstack([G, np.zeros((G.shape[0], int(m[2])), np.uint8)])
    m = re.search(r"ExtendCode \[(\d+)\] by (\d+)", body)
    if m:
        G = built[int(m[1])]
        par = np.zeros(G.shape[0], np.uint8)
        for c in range(G.shape[1]):
            par = F.add[par, G[:, c]]
        cols = [F.neg[par].reshape(-1, 1)] + [np.zeros((G.shape[0], 1), np.uint8)] * (int(m[2]) - 1)
        return np.hstack([G] + cols)
    m = re.search(r"^Dual of \[(\d+)\]", body)
    if m:
        return dual(F, built[int(m[1])])
    m = re.search(r"Puncturing of \[(\d+)\] at (\{[^}]*\})", body)
    if m:
        G = built[int(m[1])]
        pos = positions(m[2])
        return np.delete(G, pos, axis=1)
    m = re.search(r"Shortening of \[(\d+)\] at (\{[^}]*\})", body)
    if m:
        G = built[int(m[1])]
        return shorten(F, G, positions(m[2]))
    raise UnsupportedStep(desc + " | " + body[:80])


def build_entry(q, n, k):
    """Generator matrix (k x n) of the table's code for [n,k], plus its claimed d."""
    F = GF(q)
    steps = parse(fetch(q, n, k))
    if not steps:
        raise UnsupportedStep("no construction")
    built = {}
    for st in steps:
        G = build_step(F, st, built)
        G, _ = row_reduce(F, G)
        if G.shape != (st["k"], st["n"]):
            raise UnsupportedStep(f"step {st['idx']} built {G.shape}, expected {(st['k'], st['n'])}")
        built[st["idx"]] = G
    last = steps[-1]
    if (last["n"], last["k"]) != (n, k):
        raise UnsupportedStep(f"chain ends at [{last['n']},{last['k']}]")
    return built[last["idx"]], last["d"]
