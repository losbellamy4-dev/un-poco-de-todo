#!/usr/bin/env python3
"""Build the problem lists of the open-codes benchmark from the scraped tables.

classical.jsonl : one line per open entry [n,k]_q of codetables.de (lower < upper)
quantum.jsonl   : one line per open qubit entry [[n,k]] (lower < upper)
Each problem asks for a code with d >= lower + 1 and carries an estimate of the
cost of *proving* that distance (Brouwer-Zimmermann enumeration), used for tiers.
"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "codes" / "scripts"))
from qt import bz_cost, load_bounds  # noqa: E402

TIERS = [(1e8, "A"), (1e10, "B"), (1e12, "C"), (float("inf"), "D")]


def tier(cost):
    return next(t for lim, t in TIERS if cost <= lim)


def quantum_cost(n, k, d):
    """Rough B-Z cost for the normalizer (rank n+k) of an [[n,k]] code to prove distance d."""
    h = (n + k + 1) // 2
    free2 = max(0, (n + k) - 2 * (n - h)) if n - h < h else 0  # deficiency of the 2nd set
    total = 0.0
    for w in range(1, h + 1):
        total += 2 * math.comb(h, w) * 3 ** w
        lb = (w + 1) + max(0, w + 1 - k)
        if lb >= d:
            return total
    return total


def main():
    n_cl = 0
    with open(ROOT / "problems" / "classical.jsonl", "w") as f:
        for q in (2, 3, 4, 5, 7, 8, 9):
            B = load_bounds(q)
            for (n, k), (lo, hi) in sorted(B.items()):
                if lo >= hi:
                    continue
                cost = bz_cost(q, n, k, lo + 1)
                rec = dict(id=f"C{q}-{n}-{k}", q=q, n=n, k=k, record=lo, upper=hi, target=lo + 1,
                           verify_cost=float(f"{cost:.3g}"), tier=tier(cost))
                f.write(json.dumps(rec) + "\n")
                n_cl += 1
    n_q = 0
    with open(ROOT / "problems" / "quantum.jsonl", "w") as f:
        for line in open(ROOT.parent / "codes" / "data" / "qbounds_q4.tsv").readlines()[1:]:
            n, k, lo, hi = map(int, line.split())
            if lo >= hi:
                continue
            cost = quantum_cost(n, k, lo + 1)
            rec = dict(id=f"Q-{n}-{k}", n=n, k=k, record=lo, upper=hi, target=lo + 1,
                       verify_cost=float(f"{cost:.3g}"), tier=tier(cost) if n <= 64 else "D")
            f.write(json.dumps(rec) + "\n")
            n_q += 1
    print(f"classical: {n_cl} problems, quantum: {n_q} problems")


if __name__ == "__main__":
    main()
