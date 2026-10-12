import sys, time
import numpy as np
from gf import GF, min_distance, brute_min_distance

rng = np.random.default_rng(1)
fails = 0
for q in [2, 3, 4, 5, 7, 8, 9]:
    F = GF(q)
    # field axioms sanity
    for a in range(q):
        assert F.add[a, 0] == a and F.mul[a, 1] == a
        if a: assert F.mul[a, F.inv[a]] == 1
    for a in range(q):
        for b in range(q):
            for c in range(q):
                assert F.mul[a, F.add[b, c]] == F.add[F.mul[a, b], F.mul[a, c]]
                assert F.mul[a, F.mul[b, c]] == F.mul[F.mul[a, b], c]
    ntests = 0
    for trial in range(60):
        k = int(rng.integers(1, 6 if q > 4 else 8))
        if q ** k > 40000: continue
        n = int(rng.integers(k, k + 25))
        G = rng.integers(0, q, size=(k, n)).astype(np.uint8)
        if trial % 5 == 0 and k > 1:  # make rank-deficient sometimes
            G[-1] = F.add[G[0], F.mul[2 % q or 1, G[1]]]
        if trial % 7 == 0:  # sparse matrices -> low-weight words
            G[rng.random(G.shape) < 0.6] = 0
        d_ref = brute_min_distance(q, G)
        d = min_distance(q, G)
        d1 = min_distance(q, G, threads=1)
        ok = d == d_ref == d1
        # target mode consistency
        if d_ref > 0:
            t_hi = min_distance(q, G, target=d_ref + 1)
            t_eq = min_distance(q, G, target=d_ref)
            ok = ok and t_hi <= d_ref and t_eq >= d_ref
        if not ok:
            fails += 1
            print("MISMATCH q", q, "k", k, "n", n, d, d1, d_ref, t_hi, t_eq)
        ntests += 1
    print(f"q={q}: {ntests} random codes checked")
print("FAILS:", fails)
sys.exit(1 if fails else 0)
