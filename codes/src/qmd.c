/*
 * Minimum symplectic weight of additive codes over GF(4) (stabilizer formalism), n <= 64.
 *
 * A codeword is a pair (x, z) of n-bit masks; its weight is popcount(x | z)
 * (number of qubits where the Pauli is not the identity).  The code is the
 * GF(2)-span of r generators.
 *
 * Brouwer-Zimmermann over qubit positions: for an information set I_j the
 * generators are reduced so that every pivot bit lies in I_j.  The code is
 * then spanned by "units": full positions (both bits pivots: 3 nonzero
 * patterns), half positions (one pivot bit: 1 pattern) and free rows (no
 * pivot in I_j: 1 pattern).  Every codeword is a unique combination of unit
 * patterns, and its cost (number of units used) satisfies
 *   cost <= wt(word restricted to I_j) + #free rows.
 * Enumerating all combinations of cost <= w in every set proves that any
 * word not seen has weight >= sum_j max(0, w + 1 - free_j).
 *
 * Optionally a subcode S (given by generators) is excluded: words lying in S
 * are ignored (quantum distance of an impure code = min weight of N \ S).
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct { uint64_t x, z; } pv_t;

typedef struct {
    int nunits;
    int npat[192];
    pv_t pat[192][3];
    int nfree;
} qset_t;

typedef struct {
    int n, r, target;
    qset_t *sets;
    int nsets;
    volatile int stop;
    int best;
    pv_t best_word;
    long long count;
    /* subcode exclusion: membership test via reduced basis of S */
    int ns;
    pv_t sbasis[128];
    int spiv[128];          /* pivot bit index (0..2n-1) of each reduced S row */
} qctx_t;

static inline int getbit(pv_t v, int b, int n) { return b < n ? (v.x >> b) & 1 : (v.z >> (b - n)) & 1; }
static inline pv_t xorv(pv_t a, pv_t b) { pv_t c = { a.x ^ b.x, a.z ^ b.z }; return c; }
static inline int wt(pv_t v) { return __builtin_popcountll(v.x | v.z); }

static int in_subcode(const qctx_t *C, pv_t v) {
    for (int i = 0; i < C->ns; i++)
        if (getbit(v, C->spiv[i], C->n)) v = xorv(v, C->sbasis[i]);
    return v.x == 0 && v.z == 0;
}

/* Gaussian elimination: make `rows` systematic on bit columns of the given positions.
 * Returns rank; marks, per position, how many of its bits became pivots (posbits[pos]).
 * pivrow[b] = row index whose pivot is bit b, or -1. */
static int reduce_positions(pv_t *rows, int r, int n, const int *poslist, int npos, int *pivrow, int *posbits) {
    int rank = 0;
    for (int b = 0; b < 2 * n; b++) pivrow[b] = -1;
    for (int pi = 0; pi < npos && rank < r; pi++) {
        int pos = poslist[pi];
        posbits[pos] = 0;
        for (int half = 0; half < 2 && rank < r; half++) {
            int b = half ? pos + n : pos;
            int p = -1;
            for (int i = rank; i < r; i++) if (getbit(rows[i], b, n)) { p = i; break; }
            if (p < 0) continue;
            pv_t t = rows[p]; rows[p] = rows[rank]; rows[rank] = t;
            for (int i = 0; i < r; i++) if (i != rank && getbit(rows[i], b, n)) rows[i] = xorv(rows[i], rows[rank]);
            pivrow[b] = rank;
            posbits[pos]++;
            rank++;
        }
    }
    return rank;
}

/* choose exactly `left` more units from [start, nunits) and evaluate the full combinations */
static void qrec(qctx_t *C, const qset_t *S, int start, int left, pv_t acc, int *lbest, pv_t *lword, long long *cnt) {
    for (int u = start; u <= S->nunits - left && !C->stop; u++) {
        for (int p = 0; p < S->npat[u]; p++) {
            pv_t v = xorv(acc, S->pat[u][p]);
            if (left > 1) { qrec(C, S, u + 1, left - 1, v, lbest, lword, cnt); continue; }
            (*cnt)++;
            int w = wt(v);
            if (w < *lbest && w > 0) {
                if (!C->ns || !in_subcode(C, v)) {
                    *lbest = w; *lword = v;
                    if (C->target > 0 && w < C->target) C->stop = 1;
                }
            }
        }
    }
}

/* enumerate all combinations of exactly w units in set S (parallel over first unit) */
static void qenum(qctx_t *C, const qset_t *S, int w) {
    long long total = 0;
    #pragma omp parallel reduction(+:total)
    {
        int lbest = C->best;
        pv_t lword = { 0, 0 };
        long long cnt = 0;
        #pragma omp for schedule(dynamic, 1)
        for (int u0 = 0; u0 < S->nunits; u0++) {
            if (C->stop) continue;
            for (int p = 0; p < S->npat[u0]; p++) {
                pv_t v = S->pat[u0][p];
                if (w == 1) {
                    cnt++;
                    int ww = wt(v);
                    if (ww < lbest && ww > 0 && (!C->ns || !in_subcode(C, v))) {
                        lbest = ww; lword = v;
                        if (C->target > 0 && ww < C->target) C->stop = 1;
                    }
                } else {
                    qrec(C, S, u0 + 1, w - 1, v, &lbest, &lword, &cnt);
                }
            }
        }
        #pragma omp critical
        {
            if (lbest < C->best) { C->best = lbest; C->best_word = lword; }
        }
        total += cnt;
    }
    C->count += total;
}

/*
 * qmd_distance: minimum weight of nonzero words of span(G) not in span(Sub).
 *   G: r generators (X[i], Z[i]), n <= 64.  Sub: ns generators (may be 0).
 *   target > 0: early exit as soon as a word of weight < target is found
 *               (returns its weight), or returns >= target when proven.
 * Returns the distance (or bound), 0 if no qualifying nonzero word exists.
 * info[0] = words enumerated, info[1] = rank of G, info[2] = #sets, info[3] = final w.
 */
int qmd_distance(int n, int r, const uint64_t *X, const uint64_t *Z, int ns, const uint64_t *SX,
                 const uint64_t *SZ, int target, int nthreads, uint64_t *wit, long long *info) {
    if (n > 64 || r > 128 || ns > 128) return -1;
#ifdef _OPENMP
    if (nthreads > 0) omp_set_num_threads(nthreads);
#endif
    qctx_t C;
    memset(&C, 0, sizeof C);
    C.n = n; C.target = target; C.best = n + 1;
    /* reduced basis of the excluded subcode */
    if (ns) {
        pv_t sb[128];
        for (int i = 0; i < ns; i++) { sb[i].x = SX[i]; sb[i].z = SZ[i]; }
        int rk = 0;
        for (int b = 0; b < 2 * n && rk < ns; b++) {
            int p = -1;
            for (int i = rk; i < ns; i++) if (getbit(sb[i], b, n)) { p = i; break; }
            if (p < 0) continue;
            pv_t t = sb[p]; sb[p] = sb[rk]; sb[rk] = t;
            for (int i = 0; i < ns; i++) if (i != rk && getbit(sb[i], b, n)) sb[i] = xorv(sb[i], sb[rk]);
            C.spiv[rk] = b;
            rk++;
        }
        memcpy(C.sbasis, sb, sizeof(pv_t) * rk);
        C.ns = rk;
    }
    pv_t rows[128];
    for (int i = 0; i < r; i++) { rows[i].x = X[i]; rows[i].z = Z[i]; }
    /* rank of G */
    int allpos[64], posbits[64], pivrow[128];
    for (int i = 0; i < n; i++) allpos[i] = i;
    int rank = reduce_positions(rows, r, n, allpos, n, pivrow, posbits);
    r = rank;
    if (r == 0) return 0;
    qset_t *sets = calloc(64, sizeof(qset_t));
    int nsets = 0;
    char used[64] = {0};
    while (nsets < 64) {
        int pl[64], np = 0;
        for (int i = 0; i < n; i++) if (!used[i]) pl[np++] = i;
        if (!np) break;
        pv_t R[128];
        memcpy(R, rows, sizeof(pv_t) * r);
        int pb[64] = {0};
        int rk = reduce_positions(R, r, n, pl, np, pivrow, pb);
        if (rk == 0) break;
        qset_t *S = &sets[nsets++];
        S->nunits = 0;
        char isp[128] = {0};
        for (int pi = 0; pi < np; pi++) {
            int pos = pl[pi];
            if (!pb[pos]) continue;
            used[pos] = 1;
            int rx = pivrow[pos], rz = pivrow[pos + n];
            int u = S->nunits++;
            if (rx >= 0 && rz >= 0) {
                S->npat[u] = 3;
                S->pat[u][0] = R[rx]; S->pat[u][1] = R[rz]; S->pat[u][2] = xorv(R[rx], R[rz]);
                isp[rx] = isp[rz] = 1;
            } else {
                int rr = rx >= 0 ? rx : rz;
                S->npat[u] = 1; S->pat[u][0] = R[rr];
                isp[rr] = 1;
            }
        }
        S->nfree = 0;
        for (int i = 0; i < r; i++) if (!isp[i]) {
            int u = S->nunits++;
            S->npat[u] = 1; S->pat[u][0] = R[i];
            S->nfree++;
        }
        if (rk < r && S->nfree == r) { nsets--; break; }
    }
    C.sets = sets; C.nsets = nsets;
    int result = -2, w;
    int maxunits = 0;
    for (int j = 0; j < nsets; j++) if (sets[j].nunits > maxunits) maxunits = sets[j].nunits;
    int done[64] = {0};
    for (w = 1; w <= maxunits; w++) {
        for (int j = 0; j < nsets && !C.stop; j++) {
            if (w + 1 - sets[j].nfree <= 0) continue;
            for (int s = done[j] + 1; s <= w && !C.stop; s++) qenum(&C, &sets[j], s);
            done[j] = w;
        }
        if (C.stop) { result = C.best; break; }
        int lb = 0;
        for (int j = 0; j < nsets; j++) { int c = w + 1 - sets[j].nfree; if (c > 0) lb += c; }
        if (lb >= C.best) { result = C.best; break; }
        if (target > 0 && lb >= target) { result = lb; break; }
    }
    if (result == -2) result = C.best <= n ? C.best : 0;
    if (wit) { wit[0] = C.best_word.x; wit[1] = C.best_word.z; }
    if (info) { info[0] = C.count; info[1] = r; info[2] = nsets; info[3] = w; }
    free(sets);
    return result;
}
