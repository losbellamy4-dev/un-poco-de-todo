import itertools, numpy as np
from gf import GF, min_distance, min_words, extend_column
rng = np.random.default_rng(3)
nt = ne = 0
for q in [2,3,4,5,7,8,9]:
    F = GF(q)
    for trial in range(25):
        k = int(rng.integers(2, 5)); n = int(rng.integers(k+2, k+14))
        G = rng.integers(0,q,size=(k,n)).astype(np.uint8)
        if F.rank(G) < k: continue
        d = min_distance(q, G)
        # brute force list of projective weight-d words
        ref = set()
        for u in itertools.product(range(q), repeat=k):
            if not any(u): continue
            v = np.zeros(n, np.uint8)
            for c,row in zip(u,G):
                if c: v = F.add[v, F.mul[c,row]]
            if np.count_nonzero(v) == d:
                f = v[np.nonzero(v)[0][0]]
                ref.add(tuple(F.mul[F.inv[f], v]))
        W, S, piv = min_words(q, G, d)
        got = set(map(tuple, W))
        assert got == ref, (q, k, n, len(got), len(ref))
        assert min_words(q, G, d+1)[0] == -2
        # systematic matrix generates same code & messages are word[piv]
        for w in W[:5]:
            u = w[piv]
            v = np.zeros(n, np.uint8)
            for c,row in zip(u,S):
                if c: v = F.add[v, F.mul[c,row]]
            assert (v == w).all()
        nt += 1
        v = extend_column(q, W, piv, k)
        if v is not None:
            G2 = np.hstack([S, v.reshape(k,1)])
            assert min_distance(q, G2) >= d+1
            ne += 1
print("minwords tests", nt, "extensions verified", ne)
