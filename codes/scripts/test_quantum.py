import random
from quantum import *
def pauli(s):
    x = z = 0
    for i, c in enumerate(s):
        if c in "XY": x |= 1 << i
        if c in "ZY": z |= 1 << i
    return (x, z)
# [[5,1,3]]
five = [pauli(s) for s in ["XZZXI", "IXZZX", "XIXZZ", "ZXIXZ"]]
assert all(sym(a, b) == 0 for a in five for b in five)
print("[[5,1,3]] d =", quantum_distance(five, 5), brute_quantum_distance(five, 5))
steane = [pauli(s) for s in ["IIIXXXX", "IXXIIXX", "XIXIXIX", "IIIZZZZ", "IZZIIZZ", "ZIZIZIZ"]]
print("Steane d =", quantum_distance(steane, 7), brute_quantum_distance(steane, 7))
# hexacode as circulant graph? check random graphs vs brute force
rng = random.Random(1)
bad = 0
for t in range(300):
    n = rng.randint(3, 13)
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                adj[i] |= 1 << j; adj[j] |= 1 << i
    g = graph_code(adj, n)
    # random k>0 code: drop some generators of the graph state
    if t % 2:
        g = g[: rng.randint(1, n - 1)]
    d1 = quantum_distance(g, n); d2 = brute_quantum_distance(g, n)
    d3 = quantum_distance(g, n, threads=1)
    if not (d1 == d2 == d3):
        bad += 1; print("MISMATCH", n, len(g), d1, d2, d3)
    # target mode
    if d2 > 0:
        assert quantum_distance(g, n, target=d2 + 1) <= d2
        assert quantum_distance(g, n, target=d2) >= d2
print("random tests done, mismatches:", bad)
