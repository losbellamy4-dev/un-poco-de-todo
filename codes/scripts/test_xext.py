import numpy as np
from gf import GF, min_distance, words_range
from xext import aux_code, try_x, avoid
# aux codes have the claimed distance
for q in [2,3,4,5,7,8,9]:
    for r in range(1,5):
        for delta in range(1,5):
            A = aux_code(q, r, delta)
            if A is None: continue
            assert GF(q).rank(A) == r
            assert min_distance(q, A) >= delta, (q,r,delta)
print("aux codes OK")
# construction X on random small codes: check resulting distance
rng = np.random.default_rng(11)
ok = 0
for q in [3,4,5,7]:
    for t in range(30):
        k = int(rng.integers(3,6)); n = int(rng.integers(k+4,k+16))
        G = rng.integers(0,q,size=(k,n)).astype(np.uint8)
        if GF(q).rank(G) < k: continue
        d = min_distance(q,G)
        for delta in (1,2,3):
            T = d + delta
            W, wts, S, piv = words_range(q, G, d, T-1)
            for r in (1,2,3):
                A = aux_code(q, r, delta)
                if A is None or r > k: continue
                V = avoid(q, W, piv, k, r, attempts=200, inner=2000)
                if V is None: continue
                from gf import matmul
                G2 = np.hstack([S, matmul(q, V, A)])
                assert min_distance(q, G2) >= T, (q,k,n,d,delta,r)
                ok += 1
print("construction X verified:", ok)
