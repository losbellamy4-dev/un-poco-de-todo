#!/usr/bin/env python3
"""Fetch (and cache) codetables.de construction pages and parse the chain.

Cached HTML: codes/data/constructions/q{q}/{n}_{k}.txt (tag-stripped text).
parse(text) -> list of steps: {"idx": i, "n":, "k":, "d":, "desc": first line, "body": following lines}
usage: constructions.py fetch q n k [n k ...]     |   constructions.py flat q [--limit N]
"""
import re
import sys
import time
import urllib.request
from pathlib import Path

from gf import ROOT
from qt import load_bounds

CACHE = ROOT / "data" / "constructions"
STEP = re.compile(r"^\[(\d+)\]: \[(\d+), (\d+), (\d+)\] (.*)$")


def fetch(q, n, k, sleep=1.0):
    path = CACHE / f"q{q}" / f"{n}_{k}.txt"
    if path.exists():
        return path.read_text()
    url = f"https://www.codetables.de/BKLC/BKLC.php?q={q}&n={n}&k={k}"
    for i in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                html = r.read().decode("latin-1")
            break
        except Exception:
            time.sleep(2 ** (i + 1))
    else:
        raise RuntimeError(url)
    text = re.sub(r"<[^>]*>", " ", html)
    text = text.replace("&nbsp;", " ").replace("&gt;", ">").replace("&lt;", "<").replace("&amp;", "&")
    lines = [re.sub(r"[ \t]+", " ", ln).strip() for ln in text.splitlines()]
    text = "\n".join(ln for ln in lines if ln)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    time.sleep(sleep)
    return text


def parse(text):
    steps, cur = [], None
    for ln in text.splitlines():
        m = STEP.match(ln)
        if m:
            cur = dict(idx=int(m[1]), n=int(m[2]), k=int(m[3]), d=int(m[4]), desc=m[5], body=[])
            steps.append(cur)
        elif cur is not None:
            if ln.startswith(("last modified", "From Brouwer", "Lb(", "Ub(", "References")):
                cur = None
                continue
            cur["body"].append(ln)
    return steps


def flat_spots(q, double=False):
    B = load_bounds(q)
    out = []
    for (n, k), (lo, hi) in sorted(B.items()):
        nb = B.get((n + 2 if double else n + 1, k))
        if nb and nb[0] <= lo and nb[1] > lo:
            out.append((n, k))
    return out


if __name__ == "__main__":
    cmd, q = sys.argv[1], int(sys.argv[2])
    if cmd == "fetch":
        args = list(map(int, sys.argv[3:]))
        for n, k in zip(args[::2], args[1::2]):
            for s in parse(fetch(q, n, k)):
                print(s["idx"], s["n"], s["k"], s["d"], s["desc"], "|", " ".join(s["body"])[:150])
    elif cmd == "flat":
        limit = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else 10**9
        spots = flat_spots(q)[:limit]
        kinds = {}
        for i, (n, k) in enumerate(spots):
            steps = parse(fetch(q, n, k))
            base = steps[0]["desc"] if steps else "?"
            key = re.sub(r"\d+", "#", base)[:60]
            kinds[key] = kinds.get(key, 0) + 1
            if i % 50 == 0:
                print(f"{i}/{len(spots)}", file=sys.stderr, flush=True)
        for key, c in sorted(kinds.items(), key=lambda x: -x[1]):
            print(c, key)
