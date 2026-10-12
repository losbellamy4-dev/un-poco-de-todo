/*
 * Minimum distance of linear codes over GF(q), q in {2,3,4,5,7,8,9}.
 *
 * Brouwer-Zimmermann enumeration over a chain of (mostly disjoint)
 * information sets, parallelised with OpenMP.
 *
 * Field elements use the canonical index a0 + p*a1 + p^2*a2 (coefficient
 * vector w.r.t. the Conway polynomial: x^2+x+1 for GF(4), x^3+x+1 for
 * GF(8), x^2+2x+2 for GF(9)).  Inside the hot loop vectors are stored in an
 * "encoded" form chosen so that field addition is cheap and vectorisable:
 *   char 2   : identity encoding, addition = XOR
 *   GF(p)    : identity encoding, addition = min(a+b, a+b-p) on uint8
 *   GF(9)    : a0 | a1<<4, nibble-wise addition mod 3
 * The zero element is encoded as 0 in every case, so weight = #nonzero bytes.
 */
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define MAXQ 9
enum { F_CHAR2 = 0, F_PRIME = 1, F_GF9 = 2 };

typedef struct {
    int q, p, m, ftype;
    uint8_t add[MAXQ][MAXQ], mul[MAXQ][MAXQ], inv[MAXQ], neg[MAXQ];
    uint8_t enc[MAXQ], dec[256];
} field_t;

static int field_init(field_t *F, int q) {
    int p, m, poly[4] = {0};
    memset(F, 0, sizeof *F);
    switch (q) {
    case 2: p = 2; m = 1; break;
    case 3: p = 3; m = 1; break;
    case 5: p = 5; m = 1; break;
    case 7: p = 7; m = 1; break;
    case 4: p = 2; m = 2; poly[0] = 1; poly[1] = 1; break;           /* x^2+x+1 */
    case 8: p = 2; m = 3; poly[0] = 1; poly[1] = 1; poly[2] = 0; break; /* x^3+x+1 */
    case 9: p = 3; m = 2; poly[0] = 2; poly[1] = 2; break;           /* x^2+2x+2 */
    default: return -1;
    }
    F->q = q; F->p = p; F->m = m;
    F->ftype = (p == 2) ? F_CHAR2 : (m == 1 ? F_PRIME : F_GF9);
    /* coefficient vectors */
    int c[MAXQ][3];
    for (int a = 0; a < q; a++) { int t = a; for (int i = 0; i < 3; i++) { c[a][i] = t % p; t /= p; } }
    for (int a = 0; a < q; a++)
        for (int b = 0; b < q; b++) {
            int s = 0, pw = 1;
            for (int i = 0; i < m; i++) { s += ((c[a][i] + c[b][i]) % p) * pw; pw *= p; }
            F->add[a][b] = s;
            /* polynomial product reduced mod monic poly x^m + poly[m-1]x^{m-1}+...+poly[0] */
            int prod[6] = {0};
            for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) prod[i + j] += c[a][i] * c[b][j];
            for (int d = 2 * m - 2; d >= m; d--) {
                int coef = prod[d] % p;
                prod[d] = 0;
                for (int i = 0; i < m; i++) prod[d - m + i] -= coef * poly[i];
            }
            s = 0; pw = 1;
            for (int i = 0; i < m; i++) { s += (((prod[i] % p) + p) % p) * pw; pw *= p; }
            F->mul[a][b] = s;
        }
    for (int a = 0; a < q; a++) {
        for (int b = 0; b < q; b++) {
            if (F->mul[a][b] == 1) F->inv[a] = b;
            if (F->add[a][b] == 0) F->neg[a] = b;
        }
    }
    memset(F->dec, 0xff, sizeof F->dec);
    for (int a = 0; a < q; a++) {
        F->enc[a] = (F->ftype == F_GF9) ? (uint8_t)(c[a][0] | (c[a][1] << 4)) : (uint8_t)a;
        F->dec[F->enc[a]] = a;
    }
    return 0;
}

/* ---- vector kernels on encoded data ---------------------------------- */

static inline uint8_t enc_add(int ftype, int p, uint8_t a, uint8_t b) {
    if (ftype == F_CHAR2) return a ^ b;
    if (ftype == F_PRIME) { uint8_t s = a + b, t = s - p; return s < t ? s : t; }
    uint8_t s = a + b;                      /* GF(9): two trits in nibbles */
    uint8_t m = ((uint8_t)(s + 0x55) >> 3) & 0x11;
    return s - 3 * m;
}

#define DEF_KERNELS(SUF, FT)                                                   \
static void vadd_##SUF(uint8_t *restrict d, const uint8_t *restrict a,         \
                       const uint8_t *restrict b, int len, int p) {            \
    for (int t = 0; t < len; t++) d[t] = enc_add(FT, p, a[t], b[t]);           \
}                                                                              \
static int wadd_##SUF(const uint8_t *restrict a, const uint8_t *restrict b,    \
                      int len, int p) {                                        \
    int cnt = 0;                                                               \
    for (int t = 0; t < len; t++) cnt += enc_add(FT, p, a[t], b[t]) != 0;      \
    return cnt;                                                                \
}

DEF_KERNELS(c2, F_CHAR2)
DEF_KERNELS(pr, F_PRIME)
DEF_KERNELS(g9, F_GF9)

static int weight(const uint8_t *a, int len) {
    int cnt = 0;
    for (int t = 0; t < len; t++) cnt += a[t] != 0;
    return cnt;
}

/* ---- information sets ------------------------------------------------ */

typedef struct {
    int rank;
    uint8_t *S;   /* scaled rows: S[(i*(q-1) + (c-1)) * L + t] = c * row_i, encoded */
} infoset_t;

/* Reduce M (k x n, canonical, row-major) on candidate columns; pivots marked in used[]. */
static int reduce_on(const field_t *F, uint8_t *M, int k, int n, char *used) {
    int r = 0;
    for (int col = 0; col < n && r < k; col++) {
        if (used[col]) continue;
        int piv = -1;
        for (int i = r; i < k; i++) if (M[i * n + col]) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != r) for (int t = 0; t < n; t++) { uint8_t x = M[r * n + t]; M[r * n + t] = M[piv * n + t]; M[piv * n + t] = x; }
        uint8_t iv = F->inv[M[r * n + col]];
        for (int t = 0; t < n; t++) M[r * n + t] = F->mul[iv][M[r * n + t]];
        for (int i = 0; i < k; i++) {
            if (i == r || !M[i * n + col]) continue;
            uint8_t f = F->neg[M[i * n + col]];
            for (int t = 0; t < n; t++) M[i * n + t] = F->add[M[i * n + t]][F->mul[f][M[r * n + t]]];
        }
        used[col] = 1;
        r++;
    }
    return r;
}

/* ---- enumeration ------------------------------------------------------- */

typedef struct {
    const field_t *F;
    int n, k, L, q1;
    int target;            /* >0: stop as soon as a word of weight < target is found */
    volatile int stop;
    int best;              /* minimum weight found so far */
    uint8_t *best_word;    /* encoded, length L */
    long long count;
    /* collection mode (md_minwords): store every word of weight collect_d, projectively */
    int collect_d, collect_hi, maxwords, nwords, overflow, below;
    uint8_t *words;        /* nwords x n canonical, first nonzero coordinate = 1 */
    int *wts;              /* weight of each stored word (optional) */
    uint64_t *hash;        /* open addressing table, size hcap (power of 2), 0 = empty */
    int *hidx, hcap;
} ctx_t;

static void collect_word(ctx_t *C, const uint8_t *enc, int wt) {
    const field_t *F = C->F;
    const int n = C->n;
    uint8_t tmp[1024];
    int t0 = 0;
    while (t0 < n && !enc[t0]) t0++;
    uint8_t s = F->inv[F->dec[enc[t0]]];
    uint64_t h = 1469598103934665603ULL;
    for (int t = 0; t < n; t++) {
        tmp[t] = F->mul[s][F->dec[enc[t]]];
        h = (h ^ tmp[t]) * 1099511628211ULL;
    }
    if (!h) h = 1;
    #pragma omp critical(collect)
    {
        if (!C->overflow) {
            int pos = (int)(h & (C->hcap - 1)), dup = 0;
            while (C->hash[pos]) {
                if (C->hash[pos] == h && !memcmp(C->words + (size_t)C->hidx[pos] * n, tmp, n)) { dup = 1; break; }
                pos = (pos + 1) & (C->hcap - 1);
            }
            if (!dup) {
                if (C->nwords >= C->maxwords) { C->overflow = 1; C->stop = 1; }
                else {
                    memcpy(C->words + (size_t)C->nwords * n, tmp, n);
                    if (C->wts) C->wts[C->nwords] = wt;
                    C->hash[pos] = h; C->hidx[pos] = C->nwords++;
                }
            }
        }
    }
}

typedef struct {
    ctx_t *C;
    const infoset_t *I;
    int w, lbest, have;
    uint8_t *buf, *lword;
    long long total;
} worker_t;

/* buf + (chosen-1)*L holds the sum of `chosen` rows; choose the remaining ones from [start,k). */
static void rec(worker_t *W, int chosen, int start) {
    ctx_t *C = W->C;
    const int k = C->k, L = C->L, q1 = C->q1, p = C->F->p, ft = C->F->ftype, w = W->w;
    const uint8_t *cur = W->buf + (size_t)(chosen - 1) * L;
    if (C->stop) return;
    if (chosen == w - 1) {
        for (int i = start; i < k; i++) {
            const uint8_t *rr = W->I->S + (size_t)(i * q1) * L;
            for (int c = 1; c <= q1; c++) {
                const uint8_t *r2 = rr + (size_t)(c - 1) * L;
                int wt = ft == F_CHAR2 ? wadd_c2(cur, r2, L, p)
                       : ft == F_PRIME ? wadd_pr(cur, r2, L, p) : wadd_g9(cur, r2, L, p);
                if (C->collect_d && wt <= C->collect_hi && wt > 0) {
                    uint8_t *tw = W->lword + L;   /* scratch */
                    if (ft == F_CHAR2) vadd_c2(tw, cur, r2, L, p);
                    else if (ft == F_PRIME) vadd_pr(tw, cur, r2, L, p);
                    else vadd_g9(tw, cur, r2, L, p);
                    if (wt < C->collect_d) { C->below = 1; C->stop = 1; }
                    else collect_word(C, tw, wt);
                }
                if (wt < W->lbest && wt > 0) {
                    W->lbest = wt; W->have = 1;
                    if (ft == F_CHAR2) vadd_c2(W->lword, cur, r2, L, p);
                    else if (ft == F_PRIME) vadd_pr(W->lword, cur, r2, L, p);
                    else vadd_g9(W->lword, cur, r2, L, p);
                    if (C->target > 0 && wt < C->target) C->stop = 1;
                }
            }
            W->total += q1;
        }
        return;
    }
    uint8_t *nxt = W->buf + (size_t)chosen * L;
    for (int i = start; i <= k - w + chosen; i++) {
        const uint8_t *rr = W->I->S + (size_t)(i * q1) * L;
        for (int c = 1; c <= q1; c++) {
            const uint8_t *r2 = rr + (size_t)(c - 1) * L;
            if (ft == F_CHAR2) vadd_c2(nxt, cur, r2, L, p);
            else if (ft == F_PRIME) vadd_pr(nxt, cur, r2, L, p);
            else vadd_g9(nxt, cur, r2, L, p);
            rec(W, chosen + 1, i + 1);
            if (C->stop) return;
        }
    }
}

static void enumerate(ctx_t *C, const infoset_t *I, int w) {
    const int k = C->k, L = C->L, q1 = C->q1;
    long long total = 0;
    #pragma omp parallel reduction(+:total)
    {
        worker_t W = { .C = C, .I = I, .w = w, .lbest = C->best, .have = 0, .total = 0 };
        W.buf = aligned_alloc(64, (size_t)(w + 1) * L);
        W.lword = aligned_alloc(64, 2 * L);
        #pragma omp for schedule(dynamic, 1)
        for (int i0 = 0; i0 <= k - w; i0++) {
            if (C->stop) continue;
            /* first chosen row has coefficient 1 (projective normalisation) */
            memcpy(W.buf, I->S + (size_t)(i0 * q1) * L, L);
            if (w == 1) {
                int wt = weight(W.buf, L);
                W.total++;
                if (C->collect_d && wt > 0 && wt <= C->collect_hi) {
                    if (wt < C->collect_d) { C->below = 1; C->stop = 1; }
                    else collect_word(C, W.buf, wt);
                }
                if (wt > 0 && wt < W.lbest) {
                    W.lbest = wt; W.have = 1; memcpy(W.lword, W.buf, L);
                    if (C->target > 0 && wt < C->target) C->stop = 1;
                }
                continue;
            }
            rec(&W, 1, i0 + 1);
        }
        #pragma omp critical
        {
            if (W.have && W.lbest < C->best) { C->best = W.lbest; memcpy(C->best_word, W.lword, L); }
        }
        total += W.total;
        free(W.buf); free(W.lword);
    }
    C->count += total;
}

/* Build the chain of information sets; returns #sets, sets rank in *kout.
 * If sys0 != NULL it receives the generator matrix systematic on set 0 (k x n)
 * and piv0 its pivot columns (row i has its identity entry in column piv0[i]). */
static int build_sets(const field_t *F, int n, int k, const uint8_t *G, int L,
                      infoset_t *sets, int *kout, uint8_t *sys0, int *piv0) {
    const int q1 = F->q - 1;
    uint8_t *M = malloc((size_t)k * n);
    memcpy(M, G, (size_t)k * n);
    char *used = calloc(n, 1);
    int ns = 0, kk = k;
    while (ns < 256) {
        int r = reduce_on(F, M, kk, n, used);
        if (ns == 0) {
            kk = r;        /* rank of G: dependent rows are now zero, at the bottom */
            if (sys0) memcpy(sys0, M, (size_t)kk * n);
            if (piv0) { int i = 0; for (int t = 0; t < n; t++) if (used[t]) piv0[i++] = t; }
            if (kk == 0) break;
        }
        if (r == 0) break;
        infoset_t *I = &sets[ns++];
        I->rank = r;
        I->S = aligned_alloc(64, (size_t)kk * q1 * L);
        memset(I->S, 0, (size_t)kk * q1 * L);
        for (int i = 0; i < kk; i++)
            for (int c = 1; c <= q1; c++) {
                uint8_t *dst = I->S + (size_t)(i * q1 + c - 1) * L;
                for (int t = 0; t < n; t++) dst[t] = F->enc[F->mul[c][M[i * n + t]]];
            }
        if (r < kk) {
            int left = 0;
            for (int t = 0; t < n; t++) left += !used[t];
            if (!left) break;
        }
    }
    free(M); free(used);
    *kout = kk;
    return ns;
}

/*
 * md_compute: minimum distance of the code generated by G (k x n, canonical
 * elements, row-major, rows need not be independent -> rank is computed).
 *   target > 0 : return value < target  => found a nonzero word of that weight (d <= value)
 *                return value >= target => proven d >= target (value is the lower bound reached,
 *                                          or the exact d if it was determined)
 *   target == 0: return exact minimum distance.
 * witness (optional, length n): minimum-weight word found (canonical elements).
 * info (optional, 4 longs): [0]=codewords enumerated, [1]=rank, [2]=#infosets, [3]=last w.
 * Returns -1 on bad field, 0 if the code is {0}.
 */
int md_compute(int q, int n, int k, const uint8_t *G, int target, int nthreads,
               uint8_t *witness, long long *info) {
    field_t F;
    if (field_init(&F, q)) return -1;
#ifdef _OPENMP
    if (nthreads > 0) omp_set_num_threads(nthreads);
#endif
    const int L = (n + 63) & ~63, q1 = q - 1;
    infoset_t sets[256];
    int kk;
    int ns = build_sets(&F, n, k, G, L, sets, &kk, NULL, NULL);
    if (kk == 0) return 0;
    k = kk;
    ctx_t C = { .F = &F, .n = n, .k = k, .L = L, .q1 = q1, .target = target, .stop = 0,
                .best = n + 1, .count = 0 };
    C.best_word = aligned_alloc(64, L);
    memset(C.best_word, 0, L);
    int result = -2, w, done[256] = {0};
    for (w = 1; w <= k; w++) {
        for (int j = 0; j < ns && !C.stop; j++) {
            /* a set only helps the bound once w+1 > k-rank; then all sizes <= w are needed */
            if (w + 1 <= k - sets[j].rank) continue;
            for (int s = done[j] + 1; s <= w && !C.stop; s++) enumerate(&C, &sets[j], s);
            done[j] = w;
        }
        if (C.stop) { result = C.best; break; }
        int lb = 0;
        for (int j = 0; j < ns; j++) { int c = w + 1 - (k - sets[j].rank); if (c > 0) lb += c; }
        if (lb >= C.best) { result = C.best; break; }
        if (target > 0 && lb >= target) { result = lb; break; }
        if (w == k) { result = C.best; break; }
    }
    if (witness) for (int t = 0; t < n; t++) witness[t] = F.dec[C.best_word[t]];
    if (info) { info[0] = C.count; info[1] = k; info[2] = ns; info[3] = w; }
    for (int j = 0; j < ns; j++) free(sets[j].S);
    free(C.best_word);
    return result;
}

/*
 * md_minwords: list all codewords of weight d (d must be the minimum distance),
 * one per projective class (first nonzero coordinate scaled to 1).
 *   words: maxwords x n output (canonical); sys0: k x n generator systematic on
 *   columns piv0[0..k-1] (so the message of word c is c[piv0[i]]).
 * Returns the number of words, -1 if more than maxwords, -2 if a word of
 * weight < d exists, -3 on error.  G must have full rank k.
 */
int md_words_range(int q, int n, int k, const uint8_t *G, int d, int dhi, int maxwords, int nthreads,
                   uint8_t *words, int *wts, uint8_t *sys0, int *piv0) {
    field_t F;
    if (field_init(&F, q) || n > 1024) return -3;
#ifdef _OPENMP
    if (nthreads > 0) omp_set_num_threads(nthreads);
#endif
    const int L = (n + 63) & ~63, q1 = q - 1;
    infoset_t sets[256];
    int kk;
    int ns = build_sets(&F, n, k, G, L, sets, &kk, sys0, piv0);
    if (kk != k) { for (int j = 0; j < ns; j++) free(sets[j].S); return -3; }
    ctx_t C = { .F = &F, .n = n, .k = k, .L = L, .q1 = q1, .target = 0, .stop = 0,
                .best = n + 1, .count = 0, .collect_d = d, .collect_hi = dhi, .maxwords = maxwords,
                .words = words, .wts = wts };
    C.hcap = 1; while (C.hcap < 2 * maxwords + 2) C.hcap <<= 1;
    C.hash = calloc(C.hcap, sizeof(uint64_t));
    C.hidx = calloc(C.hcap, sizeof(int));
    C.best_word = aligned_alloc(64, L);
    int done[256] = {0};
    for (int w = 1; w <= k && !C.stop; w++) {
        for (int j = 0; j < ns && !C.stop; j++) {
            if (w + 1 <= k - sets[j].rank) continue;
            for (int s = done[j] + 1; s <= w && !C.stop; s++) enumerate(&C, &sets[j], s);
            done[j] = w;
        }
        int lb = 0;
        for (int j = 0; j < ns; j++) { int c = w + 1 - (k - sets[j].rank); if (c > 0) lb += c; }
        if (lb > dhi) break;
    }
    int res = C.below ? -2 : C.overflow ? -1 : C.nwords;
    for (int j = 0; j < ns; j++) free(sets[j].S);
    free(C.hash); free(C.hidx); free(C.best_word);
    return res;
}

int md_minwords(int q, int n, int k, const uint8_t *G, int d, int maxwords, int nthreads,
                uint8_t *words, uint8_t *sys0, int *piv0) {
    return md_words_range(q, n, k, G, d, d, maxwords, nthreads, words, NULL, sys0, piv0);
}

/*
 * md_extend: find v in GF(q)^k with sum_i u_i v_i != 0 for every message
 * u_i = words[j][piv[i]] (j < nw).  Appending column sys0*v^T... i.e. the
 * column v to the systematic generator raises the weight of every listed word.
 * Exhaustive over projective points if q^k <= exhaustive_limit, else random
 * sampling with `tries` attempts.  Returns 1 and fills v if found, else 0.
 */
int md_extend(int q, int n, int k, const uint8_t *words, int nw, const int *piv,
              double exhaustive_limit, long long tries, unsigned long long seed, uint8_t *v) {
    field_t F;
    if (field_init(&F, q)) return -3;
    uint8_t *U = malloc((size_t)nw * k);
    int *order = malloc(nw * sizeof(int));
    for (int j = 0; j < nw; j++) {
        order[j] = j;
        for (int i = 0; i < k; i++) U[(size_t)j * k + i] = words[(size_t)j * n + piv[i]];
    }
    double total = 1; for (int i = 0; i < k; i++) total *= q;
    int found = 0;
    uint8_t x[64];
    uint64_t rs = seed * 0x9E3779B97F4A7C15ULL + 12345;
    long long it = 0;
    int exhaustive = total <= exhaustive_limit;
    if (exhaustive) memset(x, 0, k);
    for (;;) {
        if (exhaustive) {
            /* next vector in base q (little endian), keep projective representatives only */
            int i = 0;
            while (i < k) { if (++x[i] < q) break; x[i] = 0; i++; }
            if (i == k) break;
            int top = k - 1; while (top >= 0 && !x[top]) top--;
            if (x[top] != 1) continue;
        } else {
            if (it++ >= tries) break;
            for (int i = 0; i < k; i++) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; x[i] = rs % q; }
        }
        int ok = 1;
        for (int jj = 0; jj < nw; jj++) {
            const uint8_t *u = U + (size_t)order[jj] * k;
            uint8_t s = 0;
            for (int i = 0; i < k; i++) s = F.add[s][F.mul[u[i]][x[i]]];
            if (!s) {
                ok = 0;
                if (jj) { int t = order[jj]; memmove(order + 1, order, jj * sizeof(int)); order[0] = t; }
                break;
            }
        }
        if (ok) { memcpy(v, x, k); found = 1; break; }
    }
    free(U); free(order);
    return found;
}

/* Expose field tables for the Python side (canonical encoding). */
int md_field_tables(int q, uint8_t *add, uint8_t *mul, uint8_t *inv) {
    field_t F;
    if (field_init(&F, q)) return -1;
    for (int a = 0; a < q; a++) {
        inv[a] = F.inv[a];
        for (int b = 0; b < q; b++) { add[a * q + b] = F.add[a][b]; mul[a * q + b] = F.mul[a][b]; }
    }
    return 0;
}

static inline uint8_t dot(const field_t *F, const uint8_t *u, const uint8_t *x, int k) {
    uint8_t s = 0;
    for (int i = 0; i < k; i++) s = F->add[s][F->mul[u[i]][x[i]]];
    return s;
}

/*
 * md_extend2: two new columns v1, v2 raising d by 2:
 *   every word of weight d   needs u.v1 != 0 and u.v2 != 0,
 *   every word of weight d+1 needs u.v1 != 0 or  u.v2 != 0.
 * Candidates for v1, v2 are drawn from the pool of columns good for the
 * weight-d words (exhaustive if q^k <= exhaustive_limit, else `tries` samples),
 * keeping at most `poolmax` of them.  Returns 1 and fills v1, v2 if found.
 */
int md_extend2(int q, int n, int k, const uint8_t *words, const int *wts, int nw, const int *piv, int d,
               double exhaustive_limit, long long tries, int poolmax, unsigned long long seed,
               uint8_t *v1, uint8_t *v2) {
    field_t F;
    if (field_init(&F, q) || k > 64) return -3;
    int n0 = 0, n1 = 0;
    for (int j = 0; j < nw; j++) { if (wts[j] == d) n0++; else if (wts[j] == d + 1) n1++; }
    uint8_t *U0 = malloc((size_t)(n0 + 1) * k), *U1 = malloc((size_t)(n1 + 1) * k);
    n0 = n1 = 0;
    for (int j = 0; j < nw; j++) {
        uint8_t *dst = wts[j] == d ? U0 + (size_t)n0++ * k : wts[j] == d + 1 ? U1 + (size_t)n1++ * k : NULL;
        if (dst) for (int i = 0; i < k; i++) dst[i] = words[(size_t)j * n + piv[i]];
    }
    uint8_t *pool = malloc((size_t)poolmax * k);
    int np = 0;
    double total = 1; for (int i = 0; i < k; i++) total *= q;
    int exhaustive = total <= exhaustive_limit;
    uint8_t x[64];
    uint64_t rs = seed * 0x9E3779B97F4A7C15ULL + 777;
    long long it = 0;
    int *ord0 = malloc((n0 + 1) * sizeof(int));
    for (int j = 0; j < n0; j++) ord0[j] = j;
    if (exhaustive) memset(x, 0, k);
    while (np < poolmax) {
        if (exhaustive) {
            int i = 0;
            while (i < k) { if (++x[i] < q) break; x[i] = 0; i++; }
            if (i == k) break;
            int top = k - 1; while (top >= 0 && !x[top]) top--;
            if (x[top] != 1) continue;
        } else {
            if (it++ >= tries) break;
            for (int i = 0; i < k; i++) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; x[i] = rs % q; }
        }
        int ok = 1;
        for (int jj = 0; jj < n0; jj++) {
            if (!dot(&F, U0 + (size_t)ord0[jj] * k, x, k)) {
                ok = 0;
                if (jj) { int t = ord0[jj]; memmove(ord0 + 1, ord0, jj * sizeof(int)); ord0[0] = t; }
                break;
            }
        }
        if (ok) memcpy(pool + (size_t)np++ * k, x, k);
    }
    int found = 0;
    int *R = malloc((n1 + 1) * sizeof(int));
    for (int a = 0; a < np && !found; a++) {
        const uint8_t *va = pool + (size_t)a * k;
        int nr = 0;
        for (int j = 0; j < n1; j++) if (!dot(&F, U1 + (size_t)j * k, va, k)) R[nr++] = j;
        for (int b = a; b < np && !found; b++) {
            const uint8_t *vb = pool + (size_t)b * k;
            int ok = 1;
            for (int jj = 0; jj < nr; jj++) {
                if (!dot(&F, U1 + (size_t)R[jj] * k, vb, k)) {
                    ok = 0;
                    if (jj) { int t = R[jj]; memmove(R + 1, R, jj * sizeof(int)); R[0] = t; }
                    break;
                }
            }
            if (ok) { memcpy(v1, va, k); memcpy(v2, vb, k); found = 1; }
        }
    }
    free(U0); free(U1); free(pool); free(ord0); free(R);
    return found ? 1 : (np == 0 ? -1 : 0);
}

/* Gaussian elimination pivoting in the given column order; returns rank. */
static int reduce_order(const field_t *F, uint8_t *M, int k, int n, const int *order) {
    int r = 0;
    for (int oi = 0; oi < n && r < k; oi++) {
        int col = order[oi];
        int piv = -1;
        for (int i = r; i < k; i++) if (M[i * n + col]) { piv = i; break; }
        if (piv < 0) continue;
        if (piv != r) for (int t = 0; t < n; t++) { uint8_t x = M[r * n + t]; M[r * n + t] = M[piv * n + t]; M[piv * n + t] = x; }
        uint8_t iv = F->inv[M[r * n + col]];
        for (int t = 0; t < n; t++) M[r * n + t] = F->mul[iv][M[r * n + t]];
        for (int i = 0; i < k; i++) {
            if (i == r || !M[i * n + col]) continue;
            uint8_t f = F->neg[M[i * n + col]];
            for (int t = 0; t < n; t++) M[i * n + t] = F->add[M[i * n + t]][F->mul[f][M[r * n + t]]];
        }
        r++;
    }
    return r;
}

/*
 * md_isd: randomized search for a codeword of weight < target (information set
 * decoding, Lee-Brickell style): `rounds` random information sets, all
 * combinations of up to `w` rows in each.  G must have full rank k.
 * Returns the smallest weight < target found, or n+1 if none was found
 * (which is evidence, not proof, that d >= target).
 */
int md_isd(int q, int n, int k, const uint8_t *G, int target, int rounds, int w,
           unsigned long long seed, int nthreads) {
    field_t F;
    if (field_init(&F, q)) return -1;
#ifdef _OPENMP
    if (nthreads > 0) omp_set_num_threads(nthreads);
#endif
    const int L = (n + 63) & ~63, q1 = q - 1;
    uint8_t *M = malloc((size_t)k * n);
    int *order = malloc(n * sizeof(int));
    infoset_t I;
    I.S = aligned_alloc(64, (size_t)k * q1 * L);
    memset(I.S, 0, (size_t)k * q1 * L);
    ctx_t C = { .F = &F, .n = n, .k = k, .L = L, .q1 = q1, .target = target, .stop = 0,
                .best = n + 1, .count = 0 };
    C.best_word = aligned_alloc(64, L);
    uint64_t rs = seed * 0x9E3779B97F4A7C15ULL + 99991;
    for (int t = 0; t < n; t++) order[t] = t;
    for (int r = 0; r < rounds && !C.stop; r++) {
        for (int t = n - 1; t > 0; t--) {
            rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17;
            int u = rs % (t + 1), x = order[t]; order[t] = order[u]; order[u] = x;
        }
        memcpy(M, G, (size_t)k * n);
        if (reduce_order(&F, M, k, n, order) < k) continue;
        I.rank = k;
        for (int i = 0; i < k; i++)
            for (int c = 1; c <= q1; c++) {
                uint8_t *dst = I.S + (size_t)(i * q1 + c - 1) * L;
                for (int t = 0; t < n; t++) dst[t] = F.enc[F.mul[c][M[i * n + t]]];
            }
        for (int s = 1; s <= w && s <= k && !C.stop; s++) enumerate(&C, &I, s);
    }
    int res = C.best < target ? C.best : n + 1;
    free(M); free(order); free(I.S); free(C.best_word);
    return res;
}

/*
 * md_avoid: find V (k x r, column-major in V[t*k + i]) such that u V != 0 for
 * every listed word (u_i = words[j][piv[i]]), i.e. no listed word lies in the
 * codimension-r subcode {u : uV = 0}.  Strategy: random first r-1 columns, then
 * the last column must satisfy u.v != 0 for the words still in the kernel;
 * searched exhaustively (q^k <= exhaustive_limit) or by `inner` random tries.
 * Returns 1 on success.
 */
int md_avoid(int q, int n, int k, const uint8_t *words, int nw, const int *piv, int r,
             long long attempts, long long inner, double exhaustive_limit,
             unsigned long long seed, uint8_t *V) {
    field_t F;
    if (field_init(&F, q) || k > 64 || r < 1 || r > 8) return -3;
    uint8_t *U = malloc((size_t)(nw + 1) * k);
    for (int j = 0; j < nw; j++)
        for (int i = 0; i < k; i++) U[(size_t)j * k + i] = words[(size_t)j * n + piv[i]];
    int *R = malloc((nw + 1) * sizeof(int));
    uint64_t rs = seed * 0x9E3779B97F4A7C15ULL + 4242;
    double total = 1; for (int i = 0; i < k; i++) total *= q;
    int found = 0;
    uint8_t x[64];
    for (long long a = 0; a < attempts && !found; a++) {
        /* random first r-1 columns */
        for (int t = 0; t < r - 1; t++)
            for (int i = 0; i < k; i++) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; V[t * k + i] = rs % q; }
        int nr = 0;
        for (int j = 0; j < nw; j++) {
            int zero = 1;
            for (int t = 0; t < r - 1 && zero; t++) zero = !dot(&F, U + (size_t)j * k, V + t * k, k);
            if (zero) R[nr++] = j;
        }
        /* last column */
        int exhaustive = total <= exhaustive_limit;
        long long it = 0;
        if (exhaustive) memset(x, 0, k);
        for (;;) {
            if (exhaustive) {
                int i = 0;
                while (i < k) { if (++x[i] < q) break; x[i] = 0; i++; }
                if (i == k) break;
                int top = k - 1; while (top >= 0 && !x[top]) top--;
                if (x[top] != 1) continue;
            } else {
                if (it++ >= inner) break;
                for (int i = 0; i < k; i++) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; x[i] = rs % q; }
            }
            int ok = 1;
            for (int jj = 0; jj < nr; jj++) {
                if (!dot(&F, U + (size_t)R[jj] * k, x, k)) {
                    ok = 0;
                    if (jj) { int t = R[jj]; memmove(R + 1, R, jj * sizeof(int)); R[0] = t; }
                    break;
                }
            }
            if (ok) { memcpy(V + (r - 1) * k, x, k); found = 1; break; }
        }
        if (r == 1) break;   /* nothing random to retry */
    }
    free(U); free(R);
    return found;
}
