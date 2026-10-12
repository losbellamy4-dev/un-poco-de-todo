import itertools, numpy as np
from gf import GF, min_distance, words_range, extend2
rng = np.random.default_rng(7)
nt = ok = 0
for q in [2,3,4,5,7,8,9]:
    F = GF(q)
    for trial in range(20):
        k = int(rng.integers(2, 5)); n = int(rng.integers(k+2, k+12))
        G = rng.integers(0,q,size=(k,n)).astype(np.uint8)
        if F.rank(G) < k: continue
        d = min_distance(q, G)
        ref = {}
        for u in itertools.product(range(q), repeat=k):
            if not any(u): continue
            v = np.zeros(n, np.uint8)
            for c,row in zip(u,G):
                if c: v = F.add[v, F.mul[c,row]]
            w = int(np.count_nonzero(v))
            if d <= w <= d+1:
                f = v[np.nonzero(v)[0][0]]
                ref[tuple(F.mul[F.inv[f], v])] = w
        W, wts, S, piv = words_range(q, G, d, d+1)
        got = {tuple(w): int(x) for w, x in zip(W, wts)}
        assert got == ref, (q,k,n,len(got),len(ref))
        nt += 1
        r = extend2(q, W, wts, piv, k, d)
        if r is not None:
            G2 = np.hstack([S, r[0].reshape(k,1), r[1].reshape(k,1)])
            assert min_distance(q, G2) >= d+2, "extend2 wrong"
            ok += 1
print("range tests", nt, "2-extensions verified", ok)
