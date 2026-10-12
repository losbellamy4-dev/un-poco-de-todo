"""Univariate polynomials over GF(q) (canonical element encoding from gf.py).

Polynomials are Python lists of ints, lowest degree first, without trailing zeros
(the zero polynomial is []).  Includes Cantor-Zassenhaus factorisation.
"""
import random

from gf import GF


class PolyRing:
    def __init__(self, q):
        self.F = F = GF(q)
        self.q = q
        self.add_t = F.add.tolist()
        self.mul_t = F.mul.tolist()
        self.sub_t = F.sub.tolist()
        self.inv_t = F.inv.tolist()

    @staticmethod
    def trim(a):
        while a and a[-1] == 0:
            a.pop()
        return a

    def add(self, a, b):
        if len(a) < len(b):
            a, b = b, a
        r = list(a)
        for i, c in enumerate(b):
            r[i] = self.add_t[r[i]][c]
        return self.trim(r)

    def sub(self, a, b):
        n = max(len(a), len(b))
        r = [0] * n
        for i in range(n):
            x = a[i] if i < len(a) else 0
            y = b[i] if i < len(b) else 0
            r[i] = self.sub_t[x][y]
        return self.trim(r)

    def scale(self, c, a):
        if c == 0:
            return []
        return [self.mul_t[c][x] for x in a]

    def mul(self, a, b):
        if not a or not b:
            return []
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x == 0:
                continue
            mx = self.mul_t[x]
            for j, y in enumerate(b):
                if y:
                    r[i + j] = self.add_t[r[i + j]][mx[y]]
        return self.trim(r)

    def divmod(self, a, b):
        if not b:
            raise ZeroDivisionError
        a = list(a)
        db = len(b) - 1
        ilc = self.inv_t[b[-1]]
        qt = [0] * max(0, len(a) - db)
        for i in range(len(a) - 1, db - 1, -1):
            c = a[i]
            if c == 0:
                continue
            f = self.mul_t[c][ilc]
            qt[i - db] = f
            mf = self.mul_t[f]
            for j, y in enumerate(b):
                if y:
                    a[i - db + j] = self.sub_t[a[i - db + j]][mf[y]]
        return self.trim(qt), self.trim(a[:db])

    def mod(self, a, b):
        return self.divmod(a, b)[1]

    def monic(self, a):
        if not a:
            return a
        return self.scale(self.inv_t[a[-1]], a)

    def gcd(self, a, b):
        while b:
            a, b = b, self.mod(a, b)
        return self.monic(a)

    def powmod(self, a, e, m):
        r = [1]
        a = self.mod(a, m)
        while e:
            if e & 1:
                r = self.mod(self.mul(r, a), m)
            e >>= 1
            if e:
                a = self.mod(self.mul(a, a), m)
        return r

    def deriv(self, a):
        F = self.F
        r = []
        for i in range(1, len(a)):
            c = 0
            for _ in range(i % F.p):
                c = self.add_t[c][a[i]]
            r.append(c)
        return self.trim(r)

    # ---- factorisation (input must be squarefree) ----
    def ddf(self, f):
        """Distinct-degree factorisation of a monic squarefree f -> [(g, d)]."""
        res = []
        x = [0, 1]
        h = x
        i = 0
        f = self.monic(f)
        while len(f) - 1 >= 2 * (i + 1):
            i += 1
            h = self.powmod(h, self.q, f)
            g = self.gcd(self.sub(h, x), f)
            if len(g) > 1:
                res.append((g, i))
                f, _ = self.divmod(f, g)
                h = self.mod(h, f)
        if len(f) > 1:
            res.append((f, len(f) - 1))
        return res

    def edf(self, f, d, rng):
        n = len(f) - 1
        if n == d:
            return [f]
        F = self.F
        while True:
            a = self.trim([rng.randrange(self.q) for _ in range(n)])
            if len(a) < 2:
                continue
            if F.p == 2:
                m = {2: 1, 4: 2, 8: 3}[self.q]
                b, t = list(a), list(a)
                for _ in range(m * d - 1):
                    t = self.mod(self.mul(t, t), f)
                    b = self.add(b, t)
            else:
                e = (self.q ** d - 1) // 2
                b = self.sub(self.powmod(a, e, f), [1])
            g = self.gcd(b, f)
            if 1 < len(g) < len(f):
                h, _ = self.divmod(f, g)
                return self.edf(g, d, rng) + self.edf(self.monic(h), d, rng)

    def factor_squarefree(self, f, seed=0):
        rng = random.Random(seed)
        out = []
        for g, d in self.ddf(f):
            out.extend(self.edf(self.monic(g), d, rng))
        return sorted(out, key=lambda p: (len(p), p))

    def x_n_minus(self, n, lam):
        """x^n - lam."""
        p = [0] * (n + 1)
        p[0] = self.F.neg[lam]
        p[n] = 1
        return p
