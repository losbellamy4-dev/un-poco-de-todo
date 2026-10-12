import time, numpy as np
from gf import GF, min_distance
# Binary Golay [24,12,8]: generator [I | B]
B = ["110111000101","101110001011","011100010111","111000101101","110001011011","100010110111",
     "000101101111","001011011101","010110111001","101101110001","011011100011","111111111110"]
G = np.hstack([np.eye(12, dtype=np.uint8), np.array([[int(c) for c in r] for r in B], np.uint8)])
t=time.time(); print("Golay24 d =", min_distance(2, G), f"{time.time()-t:.3f}s")
# Reed-Solomon [q-1, k, q-k] over prime fields and GF(8),GF(9) via evaluation at all nonzero elements
for q in [5,7,8,9]:
    F = GF(q)
    pts = [e for e in range(1,q)]
    for k in [2,3,4]:
        G = np.zeros((k, len(pts)), np.uint8)
        for j,x in enumerate(pts):
            v = 1
            for i in range(k):
                G[i,j] = v; v = int(F.mul[v, x])
        d = min_distance(q, G)
        assert d == q-1-k+1, (q,k,d)
print("RS codes OK")
# throughput: random [n,k] codes, exact distance with full enumeration cost
rng = np.random.default_rng(0)
for q,n,k in [(2,120,40),(7,60,12),(5,80,14),(9,50,10),(4,90,18)]:
    G = rng.integers(0,q,size=(k,n)).astype(np.uint8)
    t=time.time(); d,w,inf = min_distance(q, G, info=True); dt=time.time()-t
    print(f"q={q} [{n},{k}] d={d} enumerated={inf['enumerated']:.3e} in {dt:.2f}s -> {inf['enumerated']/dt:.3e} words/s (w={inf['w']}, sets={inf['infosets']})")
