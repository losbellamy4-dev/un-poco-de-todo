/*
 * qtx: local search over 1-generator quasi-twisted codes aimed at a record
 *      [n + j, k, T] obtained by appending j columns (j = 0, 1, 2) to a QT
 *      code [n = m*l, k, d] with d >= T - j.
 *
 * Objective (lexicographic):
 *   d < T-j : (d, -A_d)
 *   d >= T-j: (T-j, F)  with  F = j*k*ln q + sum_{w=T-j}^{T-1} A_w ln P(Bin(j,1-1/q) >= T-w)
 * F is the log of the expected number of j-column extensions reaching T when
 * the new columns are random (A_w = projective number of weight-w words).
 * When F >= ftry an exact extension search is run (md_extend / md_extend2).
 *
 * usage: qtx q l m lam k T j secs seed patience ftry g0 ... g_{l-k}
 * Output lines: HIT (d >= T already at length n), EXT (extended code),
 *               BEST d=.. F=.. (end of each restart), STATS.
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

static int q, l, m, lam, k, n, dg, T, J;
static uint8_t add[81], mul[81], inv[9], g[256];
static uint8_t *G, *words, *sys0;
static int *wts, piv[256];
static const int MAXW = 60000;
static double logP[8];

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

typedef struct { int d; double F; int nw; } obj_t;

static int better(obj_t a, obj_t b) { return a.d > b.d || (a.d == b.d && a.F > b.F + 1e-9); }
static int geq(obj_t a, obj_t b) { return a.d > b.d || (a.d == b.d && a.F >= b.F - 1e-9); }

/* Evaluate G (already built).  `floor_d`: if d < floor_d we may stop early (returns d only). */
static obj_t evaluate(int floor_d) {
    obj_t o = { 0, -1e300, 0 };
    int d = md_compute(q, n, k, G, floor_d, 1, NULL, NULL);
    if (d < floor_d) { o.d = d; return o; }
    d = md_compute(q, n, k, G, 0, 1, NULL, NULL);
    if (d >= T) { o.d = T; o.F = 1e300; return o; }   /* direct hit */
    if (d < T - J) {
        int nw = md_minwords(q, n, k, G, d, MAXW, 1, words, sys0, piv);
        o.d = d; o.F = -(double)(nw < 0 ? MAXW : nw); o.nw = nw;
        return o;
    }
    int nw = md_words_range(q, n, k, G, d, T - 1, MAXW, 1, words, wts, sys0, piv);
    o.d = T - J; o.nw = nw;
    if (nw < 0) { o.F = -1e299; return o; }
    double F = J * k * log((double)q);
    for (int i = 0; i < nw; i++) F += logP[T - wts[i]];
    o.F = F;
    return o;
}

static void print_ext(int e, int nw, const uint8_t *v1, const uint8_t *v2) {
    printf("EXT d>=%d n=%d k=%d nw=%d matrix=", T, n + e, k, nw);
    for (int r = 0; r < k; r++) {
        for (int t = 0; t < n; t++) printf("%d", sys0[r * n + t]);
        if (e >= 1) printf("%d", v1[r]);
        if (e >= 2) printf("%d", v2[r]);
        printf("%s", r < k - 1 ? "," : "");
    }
    printf("\n");
    fflush(stdout);
}

int main(int argc, char **argv) {
    if (argc < 13) { fprintf(stderr, "usage: qtx q l m lam k T j secs seed patience ftry g...\n"); return 1; }
    q = atoi(argv[1]); l = atoi(argv[2]); m = atoi(argv[3]); lam = atoi(argv[4]); k = atoi(argv[5]);
    T = atoi(argv[6]); J = atoi(argv[7]);
    double secs = atof(argv[8]);
    rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[9]) * 0xD1B54A32D192ED03ULL;
    int patience = atoi(argv[10]);
    double ftry = atof(argv[11]);
    dg = l - k;
    if (argc != 12 + dg + 1) { fprintf(stderr, "need %d coefficients of g\n", dg + 1); return 1; }
    if (J < 0 || J > 2) { fprintf(stderr, "j must be 0..2\n"); return 1; }
    if (md_field_tables(q, add, mul, inv)) return 1;
    for (int i = 0; i <= dg; i++) g[i] = atoi(argv[12 + i]);
    n = m * l;
    /* logP[r] = ln P(Bin(J, p) >= r), p = 1 - 1/q */
    {
        double p = 1.0 - 1.0 / q;
        for (int r = 0; r <= 7; r++) {
            double s = 0;
            for (int t = r; t <= J; t++) {
                double c = 1; for (int u = 0; u < t; u++) c = c * (J - u) / (u + 1);
                s += c * pow(p, t) * pow(1 - p, J - t);
            }
            logP[r] = s > 0 ? log(s) : -1e6;
        }
    }
    G = malloc((size_t)k * n);
    words = malloc((size_t)MAXW * n);
    sys0 = malloc((size_t)k * n);
    wts = malloc(MAXW * sizeof(int));
    uint8_t *f = malloc((size_t)m * k), *fb = malloc((size_t)m * k);
    long long steps = 0, restarts = 0, hits = 0, exts = 0, ext_attempts = 0;
    obj_t global = { 0, -1e300, 0 };
    double t0 = now();
    while (now() - t0 < secs) {
        restarts++;
        for (int j = 0; j < m; j++)
            for (int i = 0; i < k; i++) f[j * k + i] = j == 0 ? (i == 0) : (uint8_t)(rnd() % q);
        build(f);
        obj_t cur = evaluate(1);
        memcpy(fb, f, (size_t)m * k);
        int stale = 0, attempted = 0;
        while (stale < patience && now() - t0 < secs) {
            if (cur.F >= 1e299) {   /* direct hit at length n */
                hits++;
                build(fb);
                printf("HIT d>=%d f=", T);
                for (int j = 1; j < m; j++) {
                    printf("%s", j > 1 ? ";" : "");
                    for (int i = 0; i < k; i++) printf("%s%d", i ? "," : "", fb[j * k + i]);
                }
                printf("\n");
                fflush(stdout);
                stale = patience;
                break;
            }
            if (J > 0 && cur.d == T - J && cur.F >= ftry && !attempted) {
                attempted = 1;
                ext_attempts++;
                build(fb);
                uint8_t v1[64], v2[64];
                int nw = md_words_range(q, n, k, G, T - J, T - 1, MAXW, 1, words, wts, sys0, piv);
                int ok = 0;
                if (nw > 0 && J == 1) ok = md_extend(q, n, k, words, nw, piv, 3e7, 4000000, rnd(), v1) == 1;
                if (nw > 0 && J == 2) ok = md_extend2(q, n, k, words, wts, nw, piv, T - 2, 3e7, 2000000, 4000, rnd(), v1, v2) == 1;
                if (ok) { exts++; print_ext(J, nw, v1, v2); stale = patience; break; }
            }
            memcpy(f, fb, (size_t)m * k);
            int jb = 1 + (int)(rnd() % (m - 1)), i = (int)(rnd() % k);
            uint8_t old = f[jb * k + i], nv;
            do nv = (uint8_t)(rnd() % q); while (nv == old);
            f[jb * k + i] = nv;
            if (rnd() % 4 == 0) {
                int j2 = 1 + (int)(rnd() % (m - 1)), i2 = (int)(rnd() % k);
                f[j2 * k + i2] = (uint8_t)(rnd() % q);
            }
            build(f);
            steps++;
            obj_t nxt = evaluate(cur.d);
            if (geq(nxt, cur)) {
                int imp = better(nxt, cur);
                cur = nxt;
                memcpy(fb, f, (size_t)m * k);
                if (imp) { stale = 0; attempted = 0; } else stale++;
            } else stale++;
        }
        if (better(cur, global)) global = cur;
        if (cur.d >= T - J - 1) { printf("BEST d=%d F=%.2f nw=%d\n", cur.d, cur.F, cur.nw); fflush(stdout); }
    }
    printf("STATS steps=%lld restarts=%lld hits=%lld exts=%lld ext_attempts=%lld best_d=%d best_F=%.2f secs=%.1f\n",
           steps, restarts, hits, exts, ext_attempts, global.d, global.F, now() - t0);
    return 0;
}
