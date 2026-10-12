#!/usr/bin/env python3
"""Independent verification of candidate record codes.

For every candidate in codes/results/candidates.jsonl (or a given file):
  1. rebuild the generator matrix from its description,
  2. compute the exact minimum distance with libmindist (exact mode, all threads),
  3. recompute it with GAP/GUAVA (third-party implementation), unless --no-gap,
  4. write verified records to codes/results/records.jsonl together with the
     generator matrix and a GAP/Magma-readable description.
"""
import argparse
import json
import subprocess
import sys
import tempfile
import time

import numpy as np

from gf import GF, ROOT, min_distance
from qt import QTFamily, load_bounds


def build(rec):
    if rec["kind"] == "QT":
        fam = QTFamily(rec["q"], rec["l"], rec["lam"])
        return fam.generator(rec["g"], rec["fs"], rec["k"])
    if rec["kind"] == "matrix":
        return np.array(rec["G"], np.uint8)
    raise ValueError(rec["kind"])


def gap_field_elements(q):
    """GAP expressions for canonical elements 0..q-1 (Conway polynomial basis)."""
    F = GF(q)
    if q in (2, 3, 5, 7):
        return [f"{a}*Z({q})^0" if a else f"0*Z({q})" for a in range(q)]
    # canonical a = a0 + p*a1 (+p^2*a2): element sum a_i * w^i, w = Z(q) root of Conway poly
    p = F.p
    out = []
    for a in range(q):
        coeffs, t = [], a
        while t:
            coeffs.append(t % p)
            t //= p
        terms = [f"{c}*Z({q})^{i}" for i, c in enumerate(coeffs) if c]
        out.append("+".join(terms) if terms else f"0*Z({q})")
    return out


def gap_min_distance(q, G, timeout=3600):
    els = gap_field_elements(q)
    rows = ",\n".join("[" + ",".join(els[x] for x in row) + "]" for row in G)
    script = f"""LoadPackage("guava");
G := [{rows}];
C := GeneratorMatCode(G, GF({q}));
Print("DIM ", Dimension(C), "\\n");
Print("MINDIST ", MinimumDistance(C), "\\n");
QUIT;
"""
    with tempfile.NamedTemporaryFile("w", suffix=".g", delete=False) as f:
        f.write(script)
        path = f.name
    t = time.time()
    out = subprocess.run(["gap", "-q", "-o", "8g", path], capture_output=True, text=True, timeout=timeout).stdout
    dim = int(out.split("DIM ")[1].split()[0])
    d = int(out.split("MINDIST ")[1].split()[0])
    return dim, d, time.time() - t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("file", nargs="?", default=str(ROOT / "results" / "candidates.jsonl"))
    ap.add_argument("--no-gap", action="store_true")
    ap.add_argument("--gap-timeout", type=int, default=3600)
    args = ap.parse_args()
    out_path = ROOT / "results" / "records.jsonl"
    done = set()
    if out_path.exists():
        for line in open(out_path):
            r = json.loads(line)
            done.add((r["q"], r["n"], r["k"], r["d"]))
    seen = set()
    for line in open(args.file):
        rec = json.loads(line)
        q = rec["q"]
        claim = (q, rec["n"], rec["k"], rec.get("d_lower_bound"))
        if claim in seen or any(key[:3] == claim[:3] and key[3] >= (claim[3] or 0) for key in done):
            continue
        seen.add(claim)
        G = build(rec)
        k, n = G.shape
        B = load_bounds(q)
        lo, hi = B[(n, k)]
        t = time.time()
        d, wit, info = min_distance(q, G, info=True)
        t_c = time.time() - t
        assert info["rank"] == k, "generator matrix is not of full rank"
        assert np.count_nonzero(wit) == d
        key = (q, n, k, d)
        status = "NEW RECORD" if d > lo else "no improvement"
        print(f"q={q} [{n},{k}] exact d={d} (C engine {t_c:.1f}s); table {lo}-{hi}: {status}")
        if d <= lo or key in done:
            continue
        res = dict(rec, d=d, table=[lo, hi], verified_c=True, G=G.tolist(), witness=wit.tolist())
        if not args.no_gap:
            gdim, gd, t_g = gap_min_distance(q, G, args.gap_timeout)
            print(f"   GAP/GUAVA: dim={gdim} d={gd} ({t_g:.1f}s)")
            if gdim != k or gd != d:
                print("   !!! MISMATCH with GAP — not recorded", file=sys.stderr)
                continue
            res["verified_gap"] = True
        with open(out_path, "a") as f:
            f.write(json.dumps(res) + "\n")
        done.add(key)


if __name__ == "__main__":
    main()
