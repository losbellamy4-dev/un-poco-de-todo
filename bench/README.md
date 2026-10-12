# open-codes-bench — 113 000 open problems with exact verifiers

A benchmark for automated discovery systems (AlphaEvolve-style LLM+search
loops, classical heuristics, RL). **Every solved instance is a genuine new
mathematical result**: a linear or quantum error-correcting code beating the
best known minimum distance in Markus Grassl's tables at
[codetables.de](https://www.codetables.de), the standard reference of the field.

Why this is a good benchmark:

- **Open, not synthetic.** The problems are the gaps lower < upper bound of the
  public tables (snapshot of 2026-10-12). Nobody knows the answers.
- **Exact, cheap verification.** A submission is a generator matrix (or a
  stabilizer); its minimum distance is computed exactly by a fast
  Brouwer–Zimmermann engine (~2·10⁸ codewords/s on 4 cores) and can be
  cross-checked with GAP/GUAVA.
- **Graded difficulty.** Each problem carries the estimated cost of *proving*
  the target distance, which bounds how expensive evaluating a candidate is.
- **Hard baselines.** Standard structured searches do not solve them (see
  below), so progress means real algorithmic or mathematical insight.

## Problems

| file | entries | content |
|---|---|---|
| `problems/classical.jsonl.gz` | 106 511 | open [n,k]_q entries, q ∈ {2,3,4,5,7,8,9} |
| `problems/quantum.jsonl` | 6 678 | open qubit stabilizer entries [[n,k]], n ≤ 128 |

Each line: `id, q, n, k, record, upper, target (= record + 1), verify_cost, tier`.

| tier | proof cost (codewords enumerated) | classical | quantum |
|---|---|---|---|
| A | ≤ 10⁸ (seconds) | 13 018 | 350 |
| B | ≤ 10¹⁰ (minutes) | 7 740 | 355 |
| C | ≤ 10¹² (hours) | 9 165 | 358 |
| D | larger, or quantum n > 64 | 76 588 | 5 615 |

Rebuild from fresh tables: `codes/scripts/scrape_bounds.py`, then `bench/build_problems.py`.

## Verifying a solution

```bash
cd codes/src && make && cd ../../bench
python3 verify_submission.py examples/classical_59_8_7.json --gap
python3 verify_submission.py examples/quantum_5_1.json --gap
```

Classical format: `{"q":7,"n":59,"k":8,"G":["1000000065…", …]}`, field elements
as canonical integers a0 + p·a1 + p²·a2 for the Conway-polynomial basis (same as
Magma/GAP `Z(q)`). Quantum format: `{"n":5,"k":1,"stabilizer":["XZZXI",…]}`.
The verdict reports the exact distance, the record and `improves: true/false`.

## Baselines (session 1, ~10 CPU-hours, 4 cores)

All on Tier A/B instances; **no problem was solved**:

| method | scope | outcome |
|---|---|---|
| random 1-generator quasi-twisted codes | 410 GF(7) targets, ~2M codes; 3 848-target survey over GF(3..9) | match the record on 94 entries (up to 10 % of samples), never exceed it |
| hill-climbing on QT generators, objective (d, −A_d) | flat/double-flat spots | reaches d = record−1 with a single orbit of minimum words |
| one/two-column extension, Construction X (r ≤ 4, MDS auxiliary codes) of QT codes | 1 628 + 2 097 targets | extensions work (e.g. [58,8,39]₇→[59,8,40]₇ in 80 % of tries) but only where the table already has that value |
| Construction X applied to the table's own record codes (rebuilt from their published constructions) | 287 GF(7) entries | 0 improvements: the tables are already closed under it |
| circulant / bi-circulant graph codes [[n,0]] | n = 48, 52, 56, 64, ~6M graphs each | 0 improvements |
| cyclic additive stabilizer codes {(c, c·s)} | 367 [[n,k]] plateaus, n ≤ 64 | best codes ≤ record (e.g. [[40,8]]: 7 vs 8) |

Per-entry numbers: `baselines_session1.jsonl`.

## If you beat a record

Verify with `--gap`, then send the generator matrix to the table maintainer
(codes [at] codetables.de) as the site requests. Please cite codetables.de as the
source of the bounds.
