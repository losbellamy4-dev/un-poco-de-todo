/*
 * qtclimb: local search over 1-generator quasi-twisted codes
 *   C = < (g, g f_2, ..., g f_m) > in (F[x]/(x^l - lam))^m,  k = l - deg g
 * maximising (d, -A_d) where A_d = number of projective minimum-weight words.
 * Each step changes one coefficient of one f_j (j >= 2); moves that keep the
 * objective equal are accepted (plateau walk).  Restart after `patience`
 * non-improving steps.
 *
 * usage: qtclimb q l m lam k lo secs seed patience g0 ... g_{l-k}
 * env QT_EXT=bits as in qtsearch (1: try one-column extension of d==lo codes,
 *                                  2: try two-column extension)
 * Output:
 *   HIT d>=D f2=... f3=...          code with d >= lo+1
 *   EXT d>=D n=N k=K nw=W matrix=   extended code (rows comma separated)
 *   BEST d=D A=A                     best (d, A_d) reached (once per restart, if d >= lo-1)
 *   STATS ...
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int md_compute(int q, int n, int k, const uint8_t *G, int target, int nthreads,
               uint8_t *witness, long long *info);
int md_field_tables(int q, uint8_t *add, uint8_t *mul, uint8_t *inv);
int md_minwords(int q, int n, int k, const uint8_t *G, int d, int maxwords, int nthreads,
                uint8_t *words, uint8_t *sys0, int *piv0);
int md_extend(int q, int n, int k, const uint8_t *words, int nw, const int *piv,
              double exhaustive_limit, long long tries, unsigned long long seed, uint8_t *v);
int md_words_range(int q, int n, int k, const uint8_t *G, int d, int dhi, int maxwords, int nthreads,
                   uint8_t *words, int *wts, uint8_t *sys0, int *piv0);
int md_extend2(int q, int n, int k, const uint8_t *words, const int *wts, int nw, const int *piv, int d,
               double exhaustive_limit, long long tries, int poolmax, unsigned long long seed,
               uint8_t *v1, uint8_t *v2);

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static int q, l, m, lam, k, n, dg;
static uint8_t add[81], mul[81], inv[9], g[256];
static uint8_t *G, *words, *sys0;
static int *wts, piv[256];
static const int MAXW = 20000;

static void build(const uint8_t *f) {
    uint8_t a[512];
    for (int j = 0; j < m; j++) {
        memset(a, 0, l);
        for (int i = 0; i <= dg; i++) {
            if (!g[i]) continue;
            for (int t = 0; t < k; t++) {
                uint8_t c = f[j * k + t];
                if (!c) continue;
                int e = i + t;
                uint8_t v = mul[g[i] * q + c];
                if (e >= l) { e -= l; v = mul[v * q + lam]; }
                a[e] = add[a[e] * q + v];
            }
        }
        for (int r = 0; r < k; r++) {
            uint8_t *row = G + (size_t)r * n + (size_t)j * l;
            for (int t = 0; t < l; t++) {
                int s = t - r;
                row[t] = s >= 0 ? a[s] : mul[a[s + l] * q + lam];
            }
        }
    }
}

/* objective: d (exact unless < floor_d) and A_d (capped).  Returns d; *A set. */
static int evaluate(int floor_d, int *A) {
    int d = md_compute(q, n, k, G, floor_d, 1, NULL, NULL);
    if (d < floor_d) { *A = MAXW; return d; }
    d = md_compute(q, n, k, G, 0, 1, NULL, NULL);
    int nw = md_minwords(q, n, k, G, d, MAXW, 1, words, sys0, piv);
    *A = nw < 0 ? MAXW : nw;
    return d;
}

static void print_f(const char *tag, int d, const uint8_t *f) {
    printf("%s d>=%d", tag, d);
    for (int j = 1; j < m; j++) {
        printf(" f%d=", j + 1);
        for (int i = 0; i < k; i++) printf("%s%d", i ? "," : "", f[j * k + i]);
    }
    printf("\n");
    fflush(stdout);
}

static void print_ext(int d, int e, int nw, const uint8_t *v1, const uint8_t *v2) {
    printf("EXT d>=%d n=%d k=%d nw=%d matrix=", d, n + e, k, nw);
    for (int r = 0; r < k; r++) {
        for (int t = 0; t < n; t++) printf("%d", sys0[r * n + t]);
        printf("%d", v1[r]);
        if (e == 2) printf("%d", v2[r]);
        printf("%s", r < k - 1 ? "," : "");
    }
    printf("\n");
    fflush(stdout);
}

int main(int argc, char **argv) {
    if (argc < 11) { fprintf(stderr, "usage: qtclimb q l m lam k lo secs seed patience g...\n"); return 1; }
    q = atoi(argv[1]); l = atoi(argv[2]); m = atoi(argv[3]); lam = atoi(argv[4]); k = atoi(argv[5]);
    int lo = atoi(argv[6]);
    double secs = atof(argv[7]);
    rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[8]) * 0xD1B54A32D192ED03ULL;
    int patience = atoi(argv[9]);
    dg = l - k;
    if (argc != 10 + dg + 1) { fprintf(stderr, "need %d coefficients of g\n", dg + 1); return 1; }
    if (md_field_tables(q, add, mul, inv)) return 1;
    for (int i = 0; i <= dg; i++) g[i] = atoi(argv[10 + i]);
    n = m * l;
    int ext = getenv("QT_EXT") ? atoi(getenv("QT_EXT")) : 0;
    G = malloc((size_t)k * n);
    words = malloc((size_t)MAXW * n);
    sys0 = malloc((size_t)k * n);
    wts = malloc(MAXW * sizeof(int));
    uint8_t *f = malloc((size_t)m * k), *fb = malloc((size_t)m * k);
    double nw_limit = ((k - 1) * log((double)q) + 3 * log(10.0)) / -log(1.0 - 1.0 / q);
    long long steps = 0, restarts = 0, hits = 0, exts = 0, evals_lo = 0;
    int global_d = 0, global_A = MAXW;
    double t0 = now();
    while (now() - t0 < secs) {
        /* random start */
        restarts++;
        for (int j = 0; j < m; j++)
            for (int i = 0; i < k; i++) f[j * k + i] = j == 0 ? (i == 0) : (uint8_t)(rnd() % q);
        build(f);
        int A, d = evaluate(1, &A);
        int stale = 0;
        memcpy(fb, f, (size_t)m * k);
        int tried_ext_for_current = 0;
        while (stale < patience && now() - t0 < secs) {
            /* report / exploit current code */
            if (d >= lo + 1) {
                hits++;
                build(fb);
                print_f("HIT", d, fb);
                break;
            }
            if (d == lo && !tried_ext_for_current && ext) {
                tried_ext_for_current = 1;
                evals_lo++;
                build(fb);
                uint8_t v1[64], v2[64];
                if ((ext & 1) && A <= nw_limit) {
                    md_minwords(q, n, k, G, d, MAXW, 1, words, sys0, piv);
                    if (md_extend(q, n, k, words, A, piv, 2e7, 2000000, rnd(), v1) == 1) {
                        exts++;
                        print_ext(d + 1, 1, A, v1, v1);
                    }
                }
                if (ext & 2) {
                    int nw = md_words_range(q, n, k, G, d, d + 1, MAXW, 1, words, wts, sys0, piv);
                    if (nw > 0 && md_extend2(q, n, k, words, wts, nw, piv, d, 2e7, 1000000, 3000, rnd(), v1, v2) == 1) {
                        exts++;
                        print_ext(d + 2, 2, nw, v1, v2);
                    }
                }
            }
            /* propose a move */
            memcpy(f, fb, (size_t)m * k);
            int j = 1 + (int)(rnd() % (m - 1)), i = (int)(rnd() % k);
            uint8_t old = f[j * k + i], nv;
            do nv = (uint8_t)(rnd() % q); while (nv == old);
            f[j * k + i] = nv;
            if (rnd() % 4 == 0) {   /* occasionally a double move */
                int j2 = 1 + (int)(rnd() % (m - 1)), i2 = (int)(rnd() % k);
                f[j2 * k + i2] = (uint8_t)(rnd() % q);
            }
            build(f);
            steps++;
            int A2, d2 = evaluate(d, &A2);
            if (d2 > d || (d2 == d && A2 <= A)) {
                int better = d2 > d || A2 < A;
                d = d2; A = A2;
                memcpy(fb, f, (size_t)m * k);
                if (better) { stale = 0; tried_ext_for_current = 0; } else stale++;
            } else stale++;
        }
        if (d > global_d || (d == global_d && A < global_A)) { global_d = d; global_A = A; }
        if (d >= lo - 1) { printf("BEST d=%d A=%d\n", d, A); fflush(stdout); }
    }
    printf("STATS steps=%lld restarts=%lld hits=%lld exts=%lld at_lo=%lld best_d=%d best_A=%d nw_limit=%.0f secs=%.1f\n",
           steps, restarts, hits, exts, evals_lo, global_d, global_A, nw_limit, now() - t0);
    return 0;
}
