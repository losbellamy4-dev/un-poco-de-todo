/*
 * arcsearch: find n points of PG(k-1, q) such that every hyperplane contains
 * at most m of them  <=>  a projective [n, k, n - m]_q code.
 *
 * Tabu local search on the point set.  State: counts c[H] of chosen points on
 * every hyperplane H; cost = sum_H max(0, c[H] - m).  A move swaps a chosen
 * point p (taken from an over-full hyperplane) for an unchosen point p'
 * (best of a random sample), with a tabu tenure on re-adding removed points.
 *
 * usage: arcsearch q k n m iters seed [restarts]
 * prints "FOUND" and the k x n generator matrix (columns = points) on success.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int md_field_tables(int q, uint8_t *add, uint8_t *mul, uint8_t *inv);
int md_compute(int q, int n, int k, const uint8_t *G, int target, int nthreads,
               uint8_t *witness, long long *info);

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }

int main(int argc, char **argv) {
    if (argc < 7) { fprintf(stderr, "usage: arcsearch q k n m iters seed [restarts]\n"); return 1; }
    int q = atoi(argv[1]), k = atoi(argv[2]), n = atoi(argv[3]), m = atoi(argv[4]);
    long long iters = atoll(argv[5]);
    rs = 0x9E3779B97F4A7C15ULL ^ (uint64_t)atoll(argv[6]) * 0xD1B54A32D192ED03ULL;
    int restarts = argc > 7 ? atoi(argv[7]) : 1;
    uint8_t add[81], mul[81], inv[9];
    if (md_field_tables(q, add, mul, inv) || k > 8) return 1;
    /* enumerate normalised points (last nonzero coordinate = 1) */
    long long total = 1; for (int i = 0; i < k; i++) total *= q;
    int P = (int)((total - 1) / (q - 1));
    uint8_t *pts = malloc((size_t)P * k);
    int np = 0;
    uint8_t v[8] = {0};
    for (long long t = 1; t < total; t++) {
        long long x = t; int top = -1;
        for (int i = 0; i < k; i++) { v[i] = x % q; x /= q; if (v[i]) top = i; }
        if (v[top] != 1) continue;
        memcpy(pts + (size_t)np * k, v, k); np++;
    }
    /* incidence: point i lies on hyperplane j iff <pts_i, pts_j> = 0 */
    int per = (int)((total / q - 1) / (q - 1));       /* points per hyperplane */
    int *inc = malloc((size_t)P * per * sizeof(int));
    int *ninc = calloc(P, sizeof(int));
    for (int i = 0; i < P; i++)
        for (int j = 0; j < P; j++) {
            uint8_t s = 0;
            for (int t = 0; t < k; t++) s = add[s * q + mul[pts[i * k + t] * q + pts[j * k + t]]];
            if (!s) inc[(size_t)i * per + ninc[i]++] = j;
        }
    int *cnt = malloc(P * sizeof(int));
    int *chosen = malloc(n * sizeof(int));
    char *in = malloc(P);
    long long *tabu = malloc(P * sizeof(long long));
    int best_overall = 1 << 30;
    for (int rr = 0; rr < restarts; rr++) {
        memset(cnt, 0, P * sizeof(int)); memset(in, 0, P); memset(tabu, 0, P * sizeof(long long));
        for (int i = 0; i < n; i++) {
            int p; do p = rnd() % P; while (in[p]);
            in[p] = 1; chosen[i] = p;
            for (int t = 0; t < per; t++) cnt[inc[(size_t)p * per + t]]++;
        }
        int cost = 0;
        for (int h = 0; h < P; h++) if (cnt[h] > m) cost += cnt[h] - m;
        int best = cost;
        for (long long it = 0; it < iters && cost > 0; it++) {
            /* pick a chosen point on an over-full hyperplane: sample chosen points, prefer high excess */
            int bi = -1, bscore = -1;
            for (int s = 0; s < 8; s++) {
                int i = rnd() % n, p = chosen[i], sc = 0;
                for (int t = 0; t < per; t++) { int h = inc[(size_t)p * per + t]; if (cnt[h] > m) sc += cnt[h] - m; }
                if (sc > bscore) { bscore = sc; bi = i; }
            }
            int p = chosen[bi];
            /* gain from removing p */
            int rem = 0;
            for (int t = 0; t < per; t++) if (cnt[inc[(size_t)p * per + t]] > m) rem++;
            for (int t = 0; t < per; t++) cnt[inc[(size_t)p * per + t]]--;
            /* best replacement among a sample */
            int bp = -1, bdelta = 1 << 30;
            for (int s = 0; s < 64; s++) {
                int p2 = rnd() % P;
                if (in[p2] || p2 == p || tabu[p2] > it) continue;
                int addc = 0;
                for (int t = 0; t < per; t++) if (cnt[inc[(size_t)p2 * per + t]] >= m) addc++;
                int delta = addc - rem;
                if (delta < bdelta || (delta == bdelta && (rnd() & 1))) { bdelta = delta; bp = p2; }
            }
            if (bp < 0) { for (int t = 0; t < per; t++) cnt[inc[(size_t)p * per + t]]++; continue; }
            for (int t = 0; t < per; t++) cnt[inc[(size_t)bp * per + t]]++;
            in[p] = 0; in[bp] = 1; chosen[bi] = bp;
            tabu[p] = it + 10 + rnd() % 20;
            cost += bdelta;
            if (cost < best) best = cost;
        }
        if (best < best_overall) best_overall = best;
        if (cost == 0) {
            uint8_t *G = malloc((size_t)k * n);
            for (int i = 0; i < n; i++) for (int t = 0; t < k; t++) G[t * n + i] = pts[(size_t)chosen[i] * k + t];
            int d = md_compute(q, n, k, G, 0, 1, NULL, NULL);
            printf("FOUND q=%d n=%d k=%d d=%d matrix=", q, n, k, d);
            for (int t = 0; t < k; t++) { for (int i = 0; i < n; i++) printf("%d", G[t * n + i]); printf("%s", t < k - 1 ? "," : "\n"); }
            fflush(stdout);
            free(G);
            break;
        }
    }
    printf("STATS best_cost=%d\n", best_overall);
    return 0;
}
