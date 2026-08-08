/* ec.h -- 256-bit field arithmetic mod p = 2^256 - 2^32 - 977 and short-Weierstrass
 * point arithmetic for curves y^2 = x^3 + B (A == 0).  The addition formulas do not
 * involve B, so the same code serves the instance curve and secp256k1 itself.
 */
#ifndef EC_H
#define EC_H
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;

typedef struct { u64 v[4]; } fe;      /* fully reduced element of F_p */
typedef struct { fe x, y; int inf; } pt;    /* affine point */
typedef struct { fe X, Y, Z; int inf; } jpt; /* jacobian point */

static const u64 FP[4] = {0xFFFFFFFEFFFFFC2FULL, 0xFFFFFFFFFFFFFFFFULL,
                          0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL};
#define FCONST 0x1000003D1ULL   /* 2^32 + 977 */

static inline int fe_cmp(const fe *a, const fe *b) {
    for (int i = 3; i >= 0; i--) if (a->v[i] != b->v[i]) return a->v[i] < b->v[i] ? -1 : 1;
    return 0;
}
static inline int fe_is_zero(const fe *a) { return (a->v[0]|a->v[1]|a->v[2]|a->v[3]) == 0; }
static inline int fe_eq(const fe *a, const fe *b) { return fe_cmp(a,b) == 0; }
static inline void fe_set_u64(fe *r, u64 x) { r->v[0]=x; r->v[1]=r->v[2]=r->v[3]=0; }

static inline void fe_sub_p(fe *r) {  /* r -= p, assumes r >= p */
    u64 borrow = 0;
    for (int i = 0; i < 4; i++) {
        u128 d = (u128)r->v[i] - FP[i] - borrow;
        r->v[i] = (u64)d;
        borrow = (d >> 64) ? 1 : 0;
    }
}

static inline void fe_add(fe *r, const fe *a, const fe *b) {
    u64 c = 0;
    for (int i = 0; i < 4; i++) {
        u128 s = (u128)a->v[i] + b->v[i] + c;
        r->v[i] = (u64)s; c = (u64)(s >> 64);
    }
    if (c) {   /* r + 2^256 ; subtract p == add FCONST and drop the 2^256 */
        u64 cc = 0;
        u128 s = (u128)r->v[0] + FCONST; r->v[0] = (u64)s; cc = (u64)(s>>64);
        for (int i = 1; i < 4 && cc; i++) { s = (u128)r->v[i] + cc; r->v[i]=(u64)s; cc=(u64)(s>>64); }
    }
    if (fe_cmp(r, (const fe*)FP) >= 0) fe_sub_p(r);
}

static inline void fe_sub(fe *r, const fe *a, const fe *b) {
    u64 borrow = 0;
    for (int i = 0; i < 4; i++) {
        u128 d = (u128)a->v[i] - b->v[i] - borrow;
        r->v[i] = (u64)d; borrow = (d >> 64) ? 1 : 0;
    }
    if (borrow) { /* add p */
        u64 c = 0;
        for (int i = 0; i < 4; i++) {
            u128 s = (u128)r->v[i] + FP[i] + c;
            r->v[i] = (u64)s; c = (u64)(s>>64);
        }
    }
}

static inline void fe_neg(fe *r, const fe *a) {
    if (fe_is_zero(a)) { fe_set_u64(r,0); return; }
    fe z; memcpy(z.v, FP, sizeof z.v);
    fe_sub(r, &z, a);
}

static void fe_reduce_wide(fe *r, const u64 lo[8]) {
    u64 e[5]; u64 c = 0;
    for (int i = 0; i < 4; i++) {
        u128 v = (u128)lo[4+i] * FCONST + lo[i] + c;
        e[i] = (u64)v; c = (u64)(v >> 64);
    }
    e[4] = c;
    u128 v = (u128)e[4] * FCONST + e[0];
    r->v[0] = (u64)v; u64 cc = (u64)(v >> 64);
    for (int i = 1; i < 4; i++) { u128 s = (u128)e[i] + cc; r->v[i] = (u64)s; cc = (u64)(s>>64); }
    while (cc) {  /* fold the 2^256 overflow back in */
        u128 s = (u128)r->v[0] + (u128)cc * FCONST; r->v[0] = (u64)s; cc = (u64)(s>>64);
        for (int i = 1; i < 4 && cc; i++) { s = (u128)r->v[i] + cc; r->v[i]=(u64)s; cc=(u64)(s>>64); }
    }
    if (fe_cmp(r, (const fe*)FP) >= 0) fe_sub_p(r);
}

static inline void fe_mul(fe *r, const fe *a, const fe *b) {
    u64 lo[8]; memset(lo, 0, sizeof lo);
    for (int i = 0; i < 4; i++) {
        u64 c = 0;
        for (int j = 0; j < 4; j++) {
            u128 v = (u128)a->v[i] * b->v[j] + lo[i+j] + c;
            lo[i+j] = (u64)v; c = (u64)(v >> 64);
        }
        lo[i+4] = c;
    }
    fe_reduce_wide(r, lo);
}
static inline void fe_sqr(fe *r, const fe *a) { fe_mul(r, a, a); }

static void fe_inv(fe *r, const fe *a) {   /* a^(p-2) */
    /* p-2 = 2^256 - 2^32 - 979 */
    static const u64 E[4] = {0xFFFFFFFEFFFFFC2DULL, 0xFFFFFFFFFFFFFFFFULL,
                             0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL};
    fe base = *a, acc; fe_set_u64(&acc, 1);
    for (int i = 0; i < 256; i++) {
        if ((E[i>>6] >> (i & 63)) & 1) fe_mul(&acc, &acc, &base);
        fe_sqr(&base, &base);
    }
    *r = acc;
}

/* batch inversion: inv[i] = 1/in[i]; zero inputs are left as zero. */
static void fe_batch_inv(fe *out, const fe *in, int m, fe *scratch) {
    fe run; fe_set_u64(&run, 1);
    for (int i = 0; i < m; i++) {
        scratch[i] = run;
        if (!fe_is_zero(&in[i])) fe_mul(&run, &run, &in[i]);
    }
    fe ri; fe_inv(&ri, &run);
    for (int i = m - 1; i >= 0; i--) {
        if (fe_is_zero(&in[i])) { fe_set_u64(&out[i], 0); continue; }
        fe t; fe_mul(&t, &ri, &scratch[i]);
        fe_mul(&ri, &ri, &in[i]);
        out[i] = t;
    }
}

/* ---------------- jacobian arithmetic (a == 0) ---------------- */
static inline void j_set_inf(jpt *r) { r->inf = 1; fe_set_u64(&r->X,1); fe_set_u64(&r->Y,1); fe_set_u64(&r->Z,0); }
static inline void j_from_affine(jpt *r, const pt *P) {
    if (P->inf) { j_set_inf(r); return; }
    r->X = P->x; r->Y = P->y; fe_set_u64(&r->Z, 1); r->inf = 0;
}

static void j_dbl(jpt *r, const jpt *P) {
    if (P->inf || fe_is_zero(&P->Y)) { j_set_inf(r); return; }
    fe A_, B_, C_, D_, t1, t2;
    fe_sqr(&A_, &P->X);                 /* X^2 */
    fe_sqr(&B_, &P->Y);                 /* Y^2 */
    fe_sqr(&C_, &B_);                   /* Y^4 */
    fe_add(&t1, &P->X, &B_); fe_sqr(&t1, &t1);
    fe_sub(&t1, &t1, &A_); fe_sub(&t1, &t1, &C_);
    fe_add(&D_, &t1, &t1);              /* D = 2*((X+B)^2 - A - C) = 4XY^2 */
    fe E_; fe_add(&E_, &A_, &A_); fe_add(&E_, &E_, &A_);   /* 3X^2 */
    fe F_; fe_sqr(&F_, &E_);
    fe X3; fe_add(&t2, &D_, &D_); fe_sub(&X3, &F_, &t2);
    fe Y3; fe_sub(&t1, &D_, &X3); fe_mul(&Y3, &E_, &t1);
    fe t3; fe_add(&t3, &C_, &C_); fe_add(&t3, &t3, &t3); fe_add(&t3, &t3, &t3); /* 8C */
    fe_sub(&Y3, &Y3, &t3);
    fe Z3; fe_mul(&Z3, &P->Y, &P->Z); fe_add(&Z3, &Z3, &Z3);
    r->X = X3; r->Y = Y3; r->Z = Z3; r->inf = fe_is_zero(&Z3);
}

static void j_add(jpt *r, const jpt *P, const jpt *Q) {
    if (P->inf) { *r = *Q; return; }
    if (Q->inf) { *r = *P; return; }
    fe Z1Z1, Z2Z2, U1, U2, S1, S2, H, I, J, R_, V, t;
    fe_sqr(&Z1Z1, &P->Z); fe_sqr(&Z2Z2, &Q->Z);
    fe_mul(&U1, &P->X, &Z2Z2); fe_mul(&U2, &Q->X, &Z1Z1);
    fe_mul(&S1, &P->Y, &Q->Z); fe_mul(&S1, &S1, &Z2Z2);
    fe_mul(&S2, &Q->Y, &P->Z); fe_mul(&S2, &S2, &Z1Z1);
    if (fe_eq(&U1,&U2)) {
        if (fe_eq(&S1,&S2)) { j_dbl(r, P); return; }
        j_set_inf(r); return;
    }
    fe_sub(&H, &U2, &U1);
    fe_add(&I, &H, &H); fe_sqr(&I, &I);
    fe_mul(&J, &H, &I);
    fe_sub(&R_, &S2, &S1); fe_add(&R_, &R_, &R_);
    fe_mul(&V, &U1, &I);
    fe X3; fe_sqr(&X3, &R_); fe_sub(&X3, &X3, &J);
    fe_add(&t, &V, &V); fe_sub(&X3, &X3, &t);
    fe Y3; fe_sub(&t, &V, &X3); fe_mul(&Y3, &R_, &t);
    fe s; fe_mul(&s, &S1, &J); fe_add(&s, &s, &s);
    fe_sub(&Y3, &Y3, &s);
    fe Z3; fe_add(&Z3, &P->Z, &Q->Z); fe_sqr(&Z3, &Z3);
    fe_sub(&Z3, &Z3, &Z1Z1); fe_sub(&Z3, &Z3, &Z2Z2); fe_mul(&Z3, &Z3, &H);
    r->X = X3; r->Y = Y3; r->Z = Z3; r->inf = fe_is_zero(&Z3);
}

static void j_to_affine(pt *r, const jpt *P) {
    if (P->inf || fe_is_zero(&P->Z)) { r->inf = 1; fe_set_u64(&r->x,0); fe_set_u64(&r->y,0); return; }
    fe zi, zi2, zi3;
    fe_inv(&zi, &P->Z); fe_sqr(&zi2, &zi); fe_mul(&zi3, &zi2, &zi);
    fe_mul(&r->x, &P->X, &zi2); fe_mul(&r->y, &P->Y, &zi3); r->inf = 0;
}

/* scalar (256-bit, little-endian limbs) times affine point */
static void pt_mul(pt *r, const u64 k[4], const pt *P) {
    jpt acc; j_set_inf(&acc);
    jpt base; j_from_affine(&base, P);
    int top = -1;
    for (int i = 255; i >= 0; i--) if ((k[i>>6] >> (i&63)) & 1) { top = i; break; }
    if (top < 0) { r->inf = 1; fe_set_u64(&r->x,0); fe_set_u64(&r->y,0); return; }
    for (int i = top; i >= 0; i--) {
        j_dbl(&acc, &acc);
        if ((k[i>>6] >> (i&63)) & 1) j_add(&acc, &acc, &base);
    }
    j_to_affine(r, &acc);
}

/* affine add, denominators supplied pre-inverted (batch).  P != +-Q assumed. */
static inline void pt_add_with_inv(pt *r, const pt *P, const pt *Q, const fe *inv_dx) {
    fe lam, t, x3;
    fe_sub(&t, &Q->y, &P->y);
    fe_mul(&lam, &t, inv_dx);
    fe_sqr(&x3, &lam); fe_sub(&x3, &x3, &P->x); fe_sub(&x3, &x3, &Q->x);
    fe_sub(&t, &P->x, &x3); fe_mul(&t, &t, &lam); fe_sub(&t, &t, &P->y);
    r->x = x3; r->y = t; r->inf = 0;
}

/* full affine add (slow path, own inversion) */
static void pt_add(pt *r, const pt *P, const pt *Q) {
    jpt a, b, c; j_from_affine(&a, P); j_from_affine(&b, Q); j_add(&c, &a, &b); j_to_affine(r, &c);
}
static void pt_neg(pt *r, const pt *P) { r->inf = P->inf; r->x = P->x; fe_neg(&r->y, &P->y); }

/* ------------ hex I/O ------------- */
static int hexval(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}
static void fe_from_hex(fe *r, const char *s) {
    u64 v[4] = {0,0,0,0};
    int L = (int)strlen(s);
    for (int i = 0; i < L; i++) {
        int d = hexval(s[i]); if (d < 0) continue;
        /* shift left 4 */
        u64 carry = 0;
        for (int j = 0; j < 4; j++) { u64 nc = v[j] >> 60; v[j] = (v[j] << 4) | carry; carry = nc; }
        v[0] |= (u64)d;
    }
    memcpy(r->v, v, sizeof v);
}
static void fe_print_hex(const fe *a, char *buf) {
    sprintf(buf, "%016llx%016llx%016llx%016llx",
            (unsigned long long)a->v[3], (unsigned long long)a->v[2],
            (unsigned long long)a->v[1], (unsigned long long)a->v[0]);
}
#endif
