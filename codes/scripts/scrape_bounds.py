#!/usr/bin/env python3
"""Download the bounds tables for linear codes [n,k,d]_q from codetables.de.

Writes codes/data/bounds_q{q}.tsv with columns: n k lower upper
Usage: scrape_bounds.py q [q ...]
"""
import re
import sys
import time
import urllib.request
from pathlib import Path

MAXN = {2: 256, 3: 243, 4: 256, 5: 130, 7: 100, 8: 130, 9: 130}
CELL = re.compile(r'title="n=(\d+), k=(\d+)"[^>]*>([^<]+)</A>')
DATA = Path(__file__).resolve().parent.parent / "data"


def fetch(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read().decode("latin-1")
        except Exception as e:  # network hiccup: back off and retry
            print(f"  retry {i+1} after error: {e}", file=sys.stderr)
            time.sleep(2 ** (i + 1))
    raise RuntimeError(f"failed to fetch {url}")


def scrape(q):
    maxn = MAXN[q]
    rows = {}
    step = 32
    for n0 in range(1, maxn + 1, step):
        n1 = min(maxn, n0 + step - 1)
        url = (f"https://www.codetables.de/BKLC/Tables.php?q={q}"
               f"&n0={n0}&n1={n1}&k0=1&k1={n1}")
        html = fetch(url)
        for n, k, cell in CELL.findall(html):
            cell = cell.strip()
            if "-" in cell:
                lo, hi = cell.split("-")
            else:
                lo = hi = cell
            rows[(int(n), int(k))] = (int(lo), int(hi))
        print(f"q={q} n={n0}..{n1}: {len(rows)} entries", file=sys.stderr)
        time.sleep(1)
    DATA.mkdir(parents=True, exist_ok=True)
    out = DATA / f"bounds_q{q}.tsv"
    with open(out, "w") as f:
        f.write("n\tk\tlower\tupper\n")
        for (n, k), (lo, hi) in sorted(rows.items()):
            f.write(f"{n}\t{k}\t{lo}\t{hi}\n")
    gaps = sum(1 for lo, hi in rows.values() if lo < hi)
    print(f"q={q}: {len(rows)} entries, {gaps} open gaps -> {out}")


if __name__ == "__main__":
    for q in map(int, sys.argv[1:] or MAXN):
        scrape(q)
