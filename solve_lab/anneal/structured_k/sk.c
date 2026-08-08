/* sk.c -- structured-scalar hunting engine.
 *
 * Every mode matches on the x-coordinate only, so a hit is "up to sign"; the
 * python driver re-verifies each candidate exactly with mul(k,G)==T.
 *
 *   selftest  <px> <py> <k-hex>
 *   bsgs      <bx> <by> <tx> <ty> <babybits> <giantcount>
 *   rational  <bx> <by> <tx> <ty> <abits> <bmax>
 *   weight    <ptsfile> <maxhalf> <capbits>   subsets of size <= maxhalf per side
 *   sweight   <ptsfile> <maxhalf> <capbits>   ditto with signed digits +-2^i
 */
#include "ec.h"
#include <time.h>

static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC,&ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }

/* ------------------- open addressing hash table ------------------- */
static u64 *HK = NULL; static u64 *HV = NULL; static u64 HMASK = 0; static u64 HCNT = 0;

static void ht_init(u64 cap_bits) {
    u64 cap = 1ULL << cap_bits;
    HK = (u64*)calloc(cap, sizeof(u64));
    HV = (u64*)malloc(cap * sizeof(u64));
    if (!HK || !HV) { fprintf(stderr, "OOM allocating table 2^%llu\n", (unsigned long long)cap_bits); exit(1); }
    HMASK = cap - 1; HCNT = 0;
}
static inline u64 mix(u64 x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static inline void ht_put(u64 key, u64 val) {
    if (HCNT * 10 > (HMASK + 1) * 9) {
        fprintf(stderr, "FATAL: hash table over 90%% full (%llu / %llu) -- increase capbits\n",
                (unsigned long long)HCNT, (unsigned long long)(HMASK + 1));
        exit(2);
    }
    if (key == 0) key = 1;
    u64 i = mix(key) & HMASK;
    while (HK[i]) { if (HK[i] == key) return; i = (i + 1) & HMASK; }
    HK[i] = key; HV[i] = val; HCNT++;
}
static inline int ht_get(u64 key, u64 *out) {
    if (key == 0) key = 1;
    u64 i = mix(key) & HMASK;
    while (HK[i]) { if (HK[i] == key) { *out = HV[i]; return 1; } i = (i + 1) & HMASK; }
    return 0;
}

/* ------------------- parallel-chain walker ------------------- */
#define NCH 2048
static fe DX[NCH], DXI[NCH], SCR[NCH];
typedef void (*cb_fn)(int chain, u64 step, const pt *P, void *ctx);

static void walk(pt *cur, int nch, const pt *D, u64 L, cb_fn cb, void *ctx) {
    for (u64 s = 0; s < L; s++) {
        for (int t = 0; t < nch; t++) cb(t, s, &cur[t], ctx);
        for (int t = 0; t < nch; t++) fe_sub(&DX[t], &D->x, &cur[t].x);
        fe_batch_inv(DXI, DX, nch, SCR);
        for (int t = 0; t < nch; t++) {
            pt r;
            if (fe_is_zero(&DX[t]) || cur[t].inf) pt_add(&r, &cur[t], D);
            else pt_add_with_inv(&r, &cur[t], D, &DXI[t]);
            cur[t] = r;
        }
    }
}

static void sc_from_u64(u64 k[4], u64 v) { k[0]=v; k[1]=k[2]=k[3]=0; }

struct tabctx { u64 base; u64 stride; };
static void cb_table(int t, u64 s, const pt *P, void *ctx) {
    struct tabctx *c = (struct tabctx*)ctx;
    if (P->inf) return;
    ht_put(P->x.v[0], c->base + (u64)t * c->stride + s);
}
struct qctx { u64 base, stride; int found; u64 qi, hit; u64 nhit; };
static void cb_query(int t, u64 s, const pt *P, void *ctx) {
    struct qctx *c = (struct qctx*)ctx;
    u64 me = c->base + (u64)t*c->stride + s;
    if (P->inf) { c->nhit++; if (!c->found) { c->found = 2; c->qi = me; c->hit = 0; } return; }
    u64 v;
    if (ht_get(P->x.v[0], &v)) {
        c->nhit++;
        if (c->nhit <= 8) { printf("HIT i=%llu j=%llu\n", (unsigned long long)me, (unsigned long long)v); fflush(stdout); }
        if (!c->found) { c->found = 1; c->qi = me; c->hit = v; }
    }
}

static void read_pt(pt *P, const char *sx, const char *sy) {
    fe_from_hex(&P->x, sx); fe_from_hex(&P->y, sy); P->inf = 0;
}
static void print_pt(const char *tag, const pt *P) {
    char a[80], b[80];
    if (P->inf) { printf("%s INF\n", tag); return; }
    fe_print_hex(&P->x, a); fe_print_hex(&P->y, b);
    printf("%s %s %s\n", tag, a, b);
}

static void build_baby(const pt *Base, u64 M) {
    int nch = NCH; u64 L = M / nch;
    if (L == 0) { nch = 1; L = M; }
    static pt cur[NCH];
    u64 sc[4]; sc_from_u64(sc, L);
    pt S; pt_mul(&S, sc, Base);
    cur[0] = *Base;
    for (int t = 1; t < nch; t++) pt_add(&cur[t], &cur[t-1], &S);
    struct tabctx c = { 1, L };
    walk(cur, nch, Base, L, cb_table, &c);
}

/* =================== weight / signed-weight MITM =================== */
static pt   GEN[512];
static int  GIDX[512];      /* bit index i of the generator (+-2^i * G) */
static int  GSGN[512];
static int  NGEN;
static int  MAXHALF;
static int  PASS;
static pt   WTGT;
static int  CURIDX[8];
static u64  NPOINTS;
static int  NHITS;

/* encode a chosen combination: 9 bits per generator index, count in bits 54.. */
static u64 encode(const int *gi, int cnt) {
    u64 v = 0;
    for (int z = 0; z < cnt; z++) v |= ((u64)(gi[z] & 0x1FF)) << (9*z);
    v |= ((u64)cnt) << 54;
    return v;
}

static void wrec(const pt *cum, int depth, int lastbit) {
    int m = 0;
    static int _dummy;  (void)_dummy;
    int *bg = (int*)malloc(sizeof(int)*NGEN);
    for (int gi = 0; gi < NGEN; gi++) if (GIDX[gi] > lastbit) bg[m++] = gi;
    if (m == 0) { free(bg); return; }
    pt *res  = (pt*)malloc(sizeof(pt)*m);
    fe *bdx  = (fe*)malloc(sizeof(fe)*m);
    fe *bdxi = (fe*)malloc(sizeof(fe)*m);
    fe *bscr = (fe*)malloc(sizeof(fe)*m);

    if (cum->inf) { for (int t = 0; t < m; t++) res[t] = GEN[bg[t]]; }
    else {
        for (int t = 0; t < m; t++) fe_sub(&bdx[t], &GEN[bg[t]].x, &cum->x);
        fe_batch_inv(bdxi, bdx, m, bscr);
        for (int t = 0; t < m; t++) {
            if (fe_is_zero(&bdx[t])) pt_add(&res[t], cum, &GEN[bg[t]]);
            else pt_add_with_inv(&res[t], cum, &GEN[bg[t]], &bdxi[t]);
        }
    }
    NPOINTS += m;

    if (PASS == 0) {
        for (int t = 0; t < m; t++) {
            if (res[t].inf) continue;
            CURIDX[depth] = bg[t];
            ht_put(res[t].x.v[0], encode(CURIDX, depth+1));
        }
    } else {
        pt *qres = (pt*)malloc(sizeof(pt)*m);
        /* q = TGT - res  =  TGT + (-res) */
        for (int t = 0; t < m; t++) fe_sub(&bdx[t], &WTGT.x, &res[t].x);
        fe_batch_inv(bdxi, bdx, m, bscr);
        for (int t = 0; t < m; t++) {
            pt nr; pt_neg(&nr, &res[t]);
            if (fe_is_zero(&bdx[t])) pt_add(&qres[t], &nr, &WTGT);
            else pt_add_with_inv(&qres[t], &nr, &WTGT, &bdxi[t]);
        }
        for (int t = 0; t < m; t++) {
            CURIDX[depth] = bg[t];
            u64 qc = encode(CURIDX, depth+1);
            if (qres[t].inf) { printf("WHITEXACT q=%llu\n", (unsigned long long)qc); fflush(stdout); NHITS++; }
            else { u64 v; if (ht_get(qres[t].x.v[0], &v)) {
                printf("WHIT q=%llu t=%llu\n", (unsigned long long)qc, (unsigned long long)v); fflush(stdout); NHITS++; } }
        }
        free(qres);
    }

    if (depth + 1 < MAXHALF) {
        for (int t = 0; t < m; t++) { CURIDX[depth] = bg[t]; wrec(&res[t], depth+1, GIDX[bg[t]]); }
    }
    free(bg); free(res); free(bdx); free(bdxi); free(bscr);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: sk <mode> ...\n"); return 1; }
    const char *mode = argv[1];
    setvbuf(stdout, NULL, _IOLBF, 0);

    if (!strcmp(mode, "selftest")) {
        fe a, b, r; char buf[80];
        fe_from_hex(&a, "0123456789abcdeffedcba98765432100f1e2d3c4b5a69788796a5b4c3d2e1f00");
        fe_from_hex(&b, "fedcba98765432100123456789abcdef00112233445566778899aabbccddeeff");
        fe_mul(&r,&a,&b); fe_print_hex(&r,buf); printf("mul %s\n", buf);
        fe_add(&r,&a,&b); fe_print_hex(&r,buf); printf("add %s\n", buf);
        fe_sub(&r,&a,&b); fe_print_hex(&r,buf); printf("sub %s\n", buf);
        fe_inv(&r,&a);    fe_print_hex(&r,buf); printf("inv %s\n", buf);
        pt P; read_pt(&P, argv[2], argv[3]);
        u64 k[4]; fe kk; fe_from_hex(&kk, argv[4]); memcpy(k, kk.v, sizeof k);
        pt Q; pt_mul(&Q, k, &P); print_pt("mul_pt", &Q);
        pt R; pt_add(&R, &P, &Q); print_pt("add_pt", &R);
        /* exercise the walker: 1..8 times P via parallel chains of length 4 */
        pt cur[2]; cur[0] = P; pt S; u64 sc[4]; sc_from_u64(sc,4); pt_mul(&S,sc,&P);
        pt_add(&cur[1], &cur[0], &S);
        for (int s = 0; s < 4; s++) {
            for (int t = 0; t < 2; t++) { char tg[32]; sprintf(tg,"walk%d", t*4+s+1); print_pt(tg,&cur[t]); }
            fe dx[2], dxi[2], scr[2];
            for (int t=0;t<2;t++) fe_sub(&dx[t], &P.x, &cur[t].x);
            fe_batch_inv(dxi,dx,2,scr);
            for (int t=0;t<2;t++){ pt r2; if(fe_is_zero(&dx[t])) pt_add(&r2,&cur[t],&P); else pt_add_with_inv(&r2,&cur[t],&P,&dxi[t]); cur[t]=r2; }
        }
        return 0;
    }

    if (!strcmp(mode, "check")) {
        /* stdin: one 64-hex-digit scalar per line.  Prints the line number of any
         * k with x(k*Base) == x(Tgt), i.e. k*Base == +-Tgt. */
        pt Base, Tgt; read_pt(&Base, argv[2], argv[3]); read_pt(&Tgt, argv[4], argv[5]);
        /* nibble comb: CB[j][d] = d * 16^j * Base, j = 0..63, d = 1..15 */
        static pt CB[64][16];
        {
            pt cur = Base;                       /* 16^j * Base */
            for (int j = 0; j < 64; j++) {
                CB[j][1] = cur;
                for (int d = 2; d < 16; d++) pt_add(&CB[j][d], &CB[j][d-1], &cur);
                for (int q = 0; q < 4; q++) pt_add(&cur, &cur, &cur);
            }
        }
        char line[160]; u64 idx = 0, hits = 0; double t0 = now();
        while (fgets(line, sizeof line, stdin)) {
            fe k; fe_from_hex(&k, line);
            jpt acc; j_set_inf(&acc);
            for (int j = 0; j < 64; j++) {
                int d = (int)((k.v[j >> 4] >> (4 * (j & 15))) & 15);
                if (!d) continue;
                jpt b; j_from_affine(&b, &CB[j][d]);
                j_add(&acc, &acc, &b);
            }
            pt R; j_to_affine(&R, &acc);
            if (!R.inf && fe_eq(&R.x, &Tgt.x)) { printf("CHIT %llu\n", (unsigned long long)idx); hits++; }
            idx++;
        }
        fprintf(stderr, "[check] %llu scalars in %.1fs (%.0f/s), %llu hits\n",
                (unsigned long long)idx, now()-t0, idx/(now()-t0+1e-9), (unsigned long long)hits);
        printf("CDONE %llu\n", (unsigned long long)idx);
        return 0;
    }

    if (!strcmp(mode, "bsgs")) {
        pt Base, Tgt; read_pt(&Base, argv[2], argv[3]); read_pt(&Tgt, argv[4], argv[5]);
        int bb = atoi(argv[6]); u64 gc = strtoull(argv[7], NULL, 10);
        u64 M = 1ULL << bb;
        u64 capbits = bb + 1; if (capbits < 12) capbits = 12;
        double t0 = now();
        ht_init(capbits);
        build_baby(&Base, M);
        double t1 = now();
        fprintf(stderr, "[bsgs] baby table 2^%d = %llu entries (%llu stored) in %.1fs\n",
                bb, (unsigned long long)M, (unsigned long long)HCNT, t1-t0);
        int nch = NCH; u64 L = gc / nch; if (L == 0) { nch = 1; L = gc; }
        u64 sc[4]; sc_from_u64(sc, M);
        pt D; pt_mul(&D, sc, &Base); pt_neg(&D, &D);
        pt Dbig; sc_from_u64(sc, L); pt_mul(&Dbig, sc, &D);
        static pt cur[NCH];
        cur[0] = Tgt;
        for (int t = 1; t < nch; t++) pt_add(&cur[t], &cur[t-1], &Dbig);
        struct qctx q = {0, L, 0, 0, 0, 0};
        walk(cur, nch, &D, L, cb_query, &q);
        double t2 = now();
        fprintf(stderr, "[bsgs] %llu giant steps (total span %llu) in %.1fs\n",
                (unsigned long long)(L*(u64)nch), (unsigned long long)(L*(u64)nch*M), t2-t1);
        if (q.found == 0) printf("NOHIT M=%llu gc=%llu\n", (unsigned long long)M, (unsigned long long)(L*(u64)nch));
        else if (q.found == 2) printf("HITINF i=%llu M=%llu\n", (unsigned long long)q.qi, (unsigned long long)M);
        printf("DONE M=%llu span=%llu\n", (unsigned long long)M, (unsigned long long)(L*(u64)nch*M));
        return 0;
    }

    if (!strcmp(mode, "bsgsmulti")) {
        /* one baby table, many targets read from stdin as "<label> <x> <y>" */
        pt Base; read_pt(&Base, argv[2], argv[3]);
        int bb = atoi(argv[4]); u64 gc = strtoull(argv[5], NULL, 10);
        u64 M = 1ULL << bb;
        double t0 = now();
        ht_init(bb + 1 < 12 ? 12 : bb + 1);
        build_baby(&Base, M);
        fprintf(stderr, "[bsgsmulti] baby 2^%d in %.1fs\n", bb, now()-t0);
        int nch = NCH; u64 L = gc / nch; if (L == 0) { nch = 1; L = gc; }
        u64 sc[4]; sc_from_u64(sc, M);
        pt D; pt_mul(&D, sc, &Base); pt_neg(&D, &D);
        pt Dbig; sc_from_u64(sc, L); pt_mul(&Dbig, sc, &D);
        char lab[128], sx[100], sy[100];
        while (scanf("%127s %99s %99s", lab, sx, sy) == 3) {
            pt Tg; read_pt(&Tg, sx, sy);
            static pt cur[NCH];
            cur[0] = Tg;
            for (int t = 1; t < nch; t++) pt_add(&cur[t], &cur[t-1], &Dbig);
            struct qctx q = {0, L, 0, 0, 0, 0};
            double s0 = now();
            walk(cur, nch, &D, L, cb_query, &q);
            printf("%s %s span=%llu M=%llu %.1fs\n", q.found ? "MHIT" : "MNOHIT", lab,
                   (unsigned long long)(L*(u64)nch*M), (unsigned long long)M, now()-s0);
            if (q.found == 1) printf("  HITDATA %s i=%llu j=%llu\n", lab,
                (unsigned long long)q.qi, (unsigned long long)q.hit);
            fflush(stdout);
        }
        printf("MDONE M=%llu span=%llu\n", (unsigned long long)M, (unsigned long long)(L*(u64)nch*M));
        return 0;
    }

    if (!strcmp(mode, "rational")) {
        pt Base, Tgt; read_pt(&Base, argv[2], argv[3]); read_pt(&Tgt, argv[4], argv[5]);
        int ab = atoi(argv[6]); u64 bmax = strtoull(argv[7], NULL, 10);
        u64 M = 1ULL << ab;
        double t0 = now();
        ht_init(ab + 1 < 12 ? 12 : ab + 1);
        build_baby(&Base, M);
        double t1 = now();
        fprintf(stderr, "[rational] a-table 1..%llu in %.1fs\n", (unsigned long long)M, t1-t0);
        int nch = NCH; u64 L = bmax / nch; if (L == 0) { nch = 1; L = bmax; }
        static pt cur[NCH];
        u64 sc[4]; sc_from_u64(sc, L);
        pt S; pt_mul(&S, sc, &Tgt);
        cur[0] = Tgt;
        for (int t = 1; t < nch; t++) pt_add(&cur[t], &cur[t-1], &S);
        struct qctx q = {1, L, 0, 0, 0, 0};
        walk(cur, nch, &Tgt, L, cb_query, &q);
        double t2 = now();
        fprintf(stderr, "[rational] b-scan 1..%llu in %.1fs\n", (unsigned long long)(L*(u64)nch), t2-t1);
        if (q.found == 0) printf("NOHIT amax=%llu bmax=%llu\n", (unsigned long long)M, (unsigned long long)(L*(u64)nch));
        printf("DONE amax=%llu bmax=%llu\n", (unsigned long long)M, (unsigned long long)(L*(u64)nch));
        return 0;
    }

    if (!strcmp(mode, "weight") || !strcmp(mode, "sweight")) {
        int signed_mode = !strcmp(mode, "sweight");
        FILE *f = fopen(argv[2], "r");
        if (!f) { perror("ptsfile"); return 1; }
        static pt PTS[300]; int np = 0;
        char sx[100], sy[100], tag[32];
        while (fscanf(f, "%31s %99s %99s", tag, sx, sy) == 3) {
            if (!strcmp(tag, "T")) read_pt(&WTGT, sx, sy);
            else read_pt(&PTS[np++], sx, sy);
        }
        fclose(f);
        MAXHALF = atoi(argv[3]);
        u64 capbits = strtoull(argv[4], NULL, 10);
        fprintf(stderr, "[%s] %d chain points, maxhalf=%d, cap=2^%llu\n", mode, np, MAXHALF, (unsigned long long)capbits);
        NGEN = 0;
        for (int i = 0; i < np; i++) { GEN[NGEN]=PTS[i]; GIDX[NGEN]=i; GSGN[NGEN]=1; NGEN++; }
        if (signed_mode) for (int i = 0; i < np; i++) { pt_neg(&GEN[NGEN], &PTS[i]); GIDX[NGEN]=i; GSGN[NGEN]=-1; NGEN++; }
        ht_init(capbits);
        pt zero; zero.inf = 1; fe_set_u64(&zero.x,0); fe_set_u64(&zero.y,0);

        for (PASS = 0; PASS < 2; PASS++) {
            double t0 = now(); NPOINTS = 0;
            if (PASS == 1) {           /* the empty subset on the query side: is T itself in the table? */
                u64 v; if (ht_get(WTGT.x.v[0], &v)) { printf("WHIT q=%llu t=%llu\n", 0ULL, (unsigned long long)v); NHITS++; }
            }
            wrec(&zero, 0, -1);
            fprintf(stderr, "[%s] pass %d: %llu points, %.1fs, table=%llu\n", mode, PASS,
                    (unsigned long long)NPOINTS, now()-t0, (unsigned long long)HCNT);
        }
        printf("WDONE gens=%d maxhalf=%d hits=%d\n", NGEN, MAXHALF, NHITS);
        return 0;
    }

    fprintf(stderr, "unknown mode %s\n", mode);
    return 1;
}
