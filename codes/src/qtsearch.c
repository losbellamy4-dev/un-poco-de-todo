/*
 * qtsearch: random sampling of 1-generator quasi-twisted codes
 *   C = < (g, g f_2, ..., g f_m) > in (F[x]/(x^l - lam))^m,  dim = l - deg g = k
 * Reports the distribution of d relative to a reference value `lo` and prints
 * every code with d >= lo + 1 (as a line starting with "HIT").
 *
 * usage: qtsearch q l m lam k lo samples seed secs g0 g1 ... g_{l-k}
 *   samples<=0 means "until secs elapsed".
 * Output (stdout):
 *   HIT d>=<D> f2=<c0,c1,...> f3=...       (canonical field elements)
 *   STATS tried=<N> below=<#d<lo-2 (or <lo unless QT_HIST is set)> m2=<#d=lo-2> m1=<#lo-1> eq=<#lo> hit=<#>=lo+1>
 */
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

int main(int argc, char **argv) {
    if (argc < 11) { fprintf(stderr, "usage: qtsearch q l m lam k lo samples seed secs g...\n"); return 1; }
    int q = atoi(argv[1]), l = atoi(argv[2]), m = atoi(argv[3]), lam = atoi(argv[4]);
    int k = atoi(argv[5]), lo = atoi(argv[6]);
    long long samples = atoll(argv[7]);
    rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[8]) * 0xD1B54A32D192ED03ULL;
    double secs = atof(argv[9]);
    int dg = l - k;
    if (argc != 10 + dg + 1) { fprintf(stderr, "need %d coefficients of g\n", dg + 1); return 1; }
    uint8_t add[81], mul[81], inv[9];
    if (md_field_tables(q, add, mul, inv)) { fprintf(stderr, "bad q\n"); return 1; }
    uint8_t g[256] = {0};
    for (int i = 0; i <= dg; i++) g[i] = atoi(argv[10 + i]);
    int n = m * l;
    uint8_t *G = malloc((size_t)k * n);
    uint8_t *f = malloc((size_t)m * k), *a = malloc(l);
    long long tried = 0, cnt[5] = {0};
    double t0 = now();
    int hist = getenv("QT_HIST") != NULL;  /* histogram mode: classify d in lo-2..lo exactly */
    /* extension mode: for codes with d == lo, try to append one column raising d to lo+1 */
    /* QT_EXT bit 1: one-column extension (record at n+1), bit 2: two-column extension (record at n+2) */
    int ext = getenv("QT_EXT") ? atoi(getenv("QT_EXT")) : 0;
    const int maxwords = (ext & 2) ? 20000 : 4000;
    int *wts = (ext & 2) ? malloc(maxwords * sizeof(int)) : NULL;
    uint8_t v2[64];
    long long ext2_tried = 0, ext2_ok = 0;
    uint8_t *words = ext ? malloc((size_t)maxwords * n) : NULL, *sys0 = ext ? malloc((size_t)k * n) : NULL;
    int piv[256];
    uint8_t v[64];
    long long ext_tried = 0, ext_ok = 0, nw_sum = 0, nw_min = 1 << 30;
    double nw_limit = 0;
    /* expected number of good columns ~ q^(k-1) (1-1/q)^N: try only while >= 1e-3 */
    {
        double lnq = 0, lnr = 0;
        lnq = __builtin_log((double)q); lnr = -__builtin_log(1.0 - 1.0 / q);
        nw_limit = ((k - 1) * lnq + 3 * 2.302585) / lnr;
    }
    for (;;) {
        if (samples > 0 ? tried >= samples : now() - t0 > secs) break;
        /* f_1 = 1, f_j random of degree < k */
        for (int j = 0; j < m; j++)
            for (int i = 0; i < k; i++) f[j * k + i] = j == 0 ? (i == 0) : (uint8_t)(rnd() % q);
        for (int j = 0; j < m; j++) {
            /* a = g * f_j mod (x^l - lam) */
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
            /* rows x^r * a */
            for (int r = 0; r < k; r++) {
                uint8_t *row = G + (size_t)r * n + (size_t)j * l;
                for (int t = 0; t < l; t++) {
                    int s = t - r;
                    row[t] = s >= 0 ? a[s] : mul[a[s + l] * q + lam];
                }
            }
        }
        tried++;
        int d;
        if (!hist) {
            d = md_compute(q, n, k, G, lo, 1, NULL, NULL);
            if (d < lo) { cnt[0]++; continue; }
        } else {
            d = md_compute(q, n, k, G, lo - 2, 1, NULL, NULL);
            if (d < lo - 2) { cnt[0]++; continue; }
        }
        d = md_compute(q, n, k, G, lo + 1, 1, NULL, NULL);
        if (d >= lo + 1) {
            cnt[4]++;
            printf("HIT d>=%d", d);
            for (int j = 1; j < m; j++) {
                printf(" f%d=", j + 1);
                for (int i = 0; i < k; i++) printf("%s%d", i ? "," : "", f[j * k + i]);
            }
            printf("\n");
            fflush(stdout);
            continue;
        }
        if (!hist) {
            cnt[3]++;
            if (ext & 2) {
                int nw = md_words_range(q, n, k, G, lo, lo + 1, maxwords, 1, words, wts, sys0, piv);
                if (nw > 0) {
                    ext2_tried++;
                    if (md_extend2(q, n, k, words, wts, nw, piv, lo, 2e7, 1000000, 3000, rnd(), v, v2) == 1) {
                        ext2_ok++;
                        printf("EXT d>=%d n=%d k=%d nw=%d matrix=", lo + 2, n + 2, k, nw);
                        for (int r = 0; r < k; r++) {
                            for (int t = 0; t < n; t++) printf("%d", sys0[r * n + t]);
                            printf("%d%d%s", v[r], v2[r], r < k - 1 ? "," : "");
                        }
                        printf("\n");
                        fflush(stdout);
                    }
                }
            }
            if (ext & 1) {
                int nw = md_minwords(q, n, k, G, lo, maxwords, 1, words, sys0, piv);
                if (nw > 0) {
                    nw_sum += nw; if (nw < nw_min) nw_min = nw;
                    if (nw <= nw_limit) {
                        ext_tried++;
                        if (md_extend(q, n, k, words, nw, piv, 2e7, 2000000, rnd(), v) == 1) {
                            ext_ok++;
                            printf("EXT d>=%d n=%d k=%d nw=%d matrix=", lo + 1, n + 1, k, nw);
                            for (int r = 0; r < k; r++) {
                                for (int t = 0; t < n; t++) printf("%d", sys0[r * n + t]);
                                printf("%d%s", v[r], r < k - 1 ? "," : "");
                            }
                            printf("\n");
                            fflush(stdout);
                        }
                    }
                } else nw_sum += maxwords;
            }
            continue;
        }   /* d == lo */
        /* d in [lo-2, lo]: find out exactly where (d < lo+1 known) */
        int de = md_compute(q, n, k, G, 0, 1, NULL, NULL);
        if (de <= lo - 2) cnt[1]++;
        else if (de == lo - 1) cnt[2]++;
        else cnt[3]++;
    }
    printf("STATS tried=%lld below=%lld m2=%lld m1=%lld eq=%lld hit=%lld secs=%.1f",
           tried, cnt[0], cnt[1], cnt[2], cnt[3], cnt[4], now() - t0);
    if (ext & 1) printf(" ext_tried=%lld ext_ok=%lld nw_min=%lld nw_avg=%.1f nw_limit=%.0f", ext_tried, ext_ok,
                    nw_min, cnt[3] ? (double)nw_sum / cnt[3] : 0.0, nw_limit);
    if (ext & 2) printf(" ext2_tried=%lld ext2_ok=%lld", ext2_tried, ext2_ok);
    printf("\n");
    return 0;
}
