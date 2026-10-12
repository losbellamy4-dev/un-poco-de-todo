#!/usr/bin/env python3
"""Independent verification of quantum (stabilizer) code candidates.

Map every qubit Pauli (x, z) to the binary triple (x, z, x^z): non-identity
Paulis become weight-2 triples, so for any additive code C,
    symplectic weight(v) = Hamming weight(image(v)) / 2.
The distance of a [[n,0]] code (= min weight of S) is therefore half the
minimum distance of the binary [3n, n] image code, which is computed with
  (1) qmd (symplectic B-Z),  (2) libmindist md_compute on the binary image,
  (3) GAP/GUAVA on the binary image.
For k > 0, (2) and (3) give d(N(S)) <= d; equality with (1) certifies a pure code.

usage: qverify.py [file]   (default results/qcandidates.jsonl; writes results/qrecords.jsonl)
"""
import json
import sys
import time

import numpy as np

from gf import ROOT, min_distance
from qcyclic import CyclicFamily
from qsearch import load_qbounds
from quantum import graph_code, normalizer, qdist, rank
from verify import gap_min_distance


def binary_image(gens, n):
    G = np.zeros((len(gens), 3 * n), np.uint8)
    for r, (x, z) in enumerate(gens):
        for i in range(n):
            a, b = (x >> i) & 1, (z >> i) & 1
            G[r, 3 * i], G[r, 3 * i + 1], G[r, 3 * i + 2] = a, b, a ^ b
    return G


def stabilizer_of(rec):
    n = rec["n"]
    if "adj" in rec:
        return graph_code(rec["adj"], n)
    if rec.get("family") == "cyclic-additive":
        return CyclicFamily(n).stabilizer(rec["g"], rec["s"])
    if "stab" in rec:
        return [tuple(v) for v in rec["stab"]]
    raise ValueError("unknown candidate format")


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "results" / "qcandidates.jsonl")
    B = load_qbounds()
    out = ROOT / "results" / "qrecords.jsonl"
    seen = set()
    for line in open(path):
        rec = json.loads(line)
        n, k = rec["n"], rec["k"]
        key = (n, k, rec["d"])
        if key in seen:
            continue
        seen.add(key)
        S = stabilizer_of(rec)
        assert rank(S, n) == n - k
        # commutation
        for a in S:
            for b in S:
                assert (bin(a[0] & b[1]).count("1") + bin(a[1] & b[0]).count("1")) % 2 == 0
        t = time.time()
        if k == 0:
            d1 = qdist(S, n)
            Gb = binary_image(S, n)
        else:
            N = normalizer(S, n)
            d1 = qdist(N, n, exclude=S)
            Gb = binary_image(N, n)
        t1 = time.time() - t
        t = time.time()
        dH = min_distance(2, Gb)
        t2 = time.time() - t
        lo, hi = B[(n, k)]
        print(f"[[{n},{k}]]: qmd d={d1} ({t1:.1f}s); binary image d_H={dH} -> {dH/2} ({t2:.1f}s); table {lo}-{hi}",
              flush=True)
        ok = (dH == 2 * d1) if k == 0 else (dH <= 2 * d1)
        res = dict(rec, d=d1, d_binary_image=dH, pure=(dH == 2 * d1), table=[lo, hi])
        if d1 > lo and ok:
            gd = gap_min_distance(2, Gb)
            print(f"   GAP/GUAVA binary image: dim={gd[0]} d={gd[1]} ({gd[2]:.1f}s)", flush=True)
            res["d_gap_binary_image"] = gd[1]
            if k == 0 and gd[1] != 2 * d1:
                print("   !!! GAP mismatch", flush=True)
                continue
            res["verified"] = True
            res["stab"] = [list(v) for v in S]
            with open(out, "a") as f:
                f.write(json.dumps(res) + "\n")
            print(f"   VERIFIED NEW RECORD [[{n},{k},{d1}]] (table lower bound {lo})", flush=True)


if __name__ == "__main__":
    main()
