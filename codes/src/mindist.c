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
} ctx_t;

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
        W.lword = aligned_alloc(64, L);
        #pragma omp for schedule(dynamic, 1)
        for (int i0 = 0; i0 <= k - w; i0++) {
            if (C->stop) continue;
            /* first chosen row has coefficient 1 (projective normalisation) */
            memcpy(W.buf, I->S + (size_t)(i0 * q1) * L, L);
            if (w == 1) {
                int wt = weight(W.buf, L);
                W.total++;
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
    uint8_t *M = malloc((size_t)k * n);
    memcpy(M, G, (size_t)k * n);
    char *used = calloc(n, 1);
    /* first pass also determines the true dimension */
    infoset_t sets[256];
    int ns = 0, kk = k;
    while (ns < 256) {
        int r = reduce_on(&F, M, kk, n, used);
        if (ns == 0) {
            kk = r;        /* rank of G: discard dependent rows (they are zero now, at bottom) */
            if (kk == 0) { free(M); free(used); return 0; }
        }
        if (r == 0) break;
        infoset_t *I = &sets[ns++];
        I->rank = r;
        I->S = aligned_alloc(64, (size_t)kk * q1 * L);
        memset(I->S, 0, (size_t)kk * q1 * L);
        for (int i = 0; i < kk; i++)
            for (int c = 1; c <= q1; c++) {
                uint8_t *dst = I->S + (size_t)(i * q1 + c - 1) * L;
                for (int t = 0; t < n; t++) dst[t] = F.enc[F.mul[c][M[i * n + t]]];
            }
        if (r < kk) {
            /* partial set: continue with remaining columns until exhausted */
            int left = 0;
            for (int t = 0; t < n; t++) left += !used[t];
            if (!left) break;
        }
    }
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
    free(C.best_word); free(M); free(used);
    return result;
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
