#!/usr/bin/env python3
"""Encode 'find selector bits b with fold(b) = T' for a ladder instance, into SMT-LIB.

Shape (identical to the reduced problem, only the prime is scaled):
  acc_0 = G  (bit 0 forced on)
  for i in 1..n-1:  S_i = acc_{i-1} + L_i   (chord law, division-free)
                    acc_i = acc_{i-1} + b_i*(S_i - acc_{i-1})
  acc_{n-1} = T
Per stage: 7 modular multiplications, 1 boolean, plus a non-degeneracy guard.
"""
import sys, json

def smt_int(d, path):
    p, n, lad, T, G = d['p'], d['n'], d['lad'], d['T'], d['G']
    L = []
    A = ['(declare-fun %s () Int)' % v for v in []]
    o = []
    w = o.append
    w('(set-logic QF_NIA)')
    decls, asserts = [], []
    def V(name):
        decls.append('(declare-fun %s () Int)' % name)
        asserts.append('(assert (and (<= 0 %s) (< %s %d)))' % (name, name, p))
        return name
    def mulmod(a, b, name):
        r = V(name); q = V(name + '_q')
        asserts.append('(assert (= (* %s %s) (+ (* %d %s) %s)))' % (a, b, p, q, r))
        return r
    def red(expr, name):
        """reduce an affine expression (true value in (-4p,4p)) mod p into a var in [0,p)"""
        r = V(name); q = name + '_r'
        decls.append('(declare-fun %s () Int)' % q)
        asserts.append('(assert (and (<= (- 4) %s) (<= %s 4)))' % (q, q))
        asserts.append('(assert (= %s (+ (* %d %s) %s)))' % (expr, p, q, r))
        return r
    ax, ay = str(G[0]), str(G[1])
    for i in range(1, n):
        bx, by = lad[i]
        bi = 'b%d' % i
        decls.append('(declare-fun %s () Int)' % bi)
        asserts.append('(assert (or (= %s 0) (= %s 1)))' % (bi, bi))
        dd = red('(- %d %s)' % (bx, ax), 'd%d' % i)
        nn = red('(- %d %s)' % (by, ay), 'n%d' % i)
        asserts.append('(assert (not (= %s 0)))' % dd)
        d2 = mulmod(dd, dd, 'd2_%d' % i)
        n2 = mulmod(nn, nn, 'n2_%d' % i)
        sx = V('sx%d' % i); sy = V('sy%d' % i)
        t = red('(+ %s %s %d)' % (sx, ax, bx), 't%d' % i)
        lhs = mulmod(t, d2, 'td%d' % i)
        asserts.append('(assert (= %s %s))' % (lhs, n2))
        u = red('(+ %s %s)' % (sy, ay), 'u%d' % i)
        v = red('(- %s %s)' % (ax, sx), 'v%d' % i)
        l1 = mulmod(u, dd, 'ud%d' % i); l2 = mulmod(nn, v, 'nv%d' % i)
        asserts.append('(assert (= %s %s))' % (l1, l2))
        # mux
        nax = V('ax%d' % i); nay = V('ay%d' % i)
        gx = mulmod(bi, red('(- %s %s)' % (sx, ax), 'dx%d' % i), 'gx%d' % i)
        gy = mulmod(bi, red('(- %s %s)' % (sy, ay), 'dy%d' % i), 'gy%d' % i)
        asserts.append('(assert (= %s %s))' % (nax, red('(+ %s %s)' % (ax, gx), 'mx%d' % i)))
        asserts.append('(assert (= %s %s))' % (nay, red('(+ %s %s)' % (ay, gy), 'my%d' % i)))
        ax, ay = nax, nay
    asserts.append('(assert (= %s %d))' % (ax, T[0]))
    asserts.append('(assert (= %s %d))' % (ay, T[1]))
    o += [x for x in decls if x] + asserts
    w('(check-sat)')
    w('(get-model)')
    open(path, 'w').write('\n'.join(o) + '\n')
    return path


def smt_bv(d, path):
    p, n, lad, T, G = d['p'], d['n'], d['lad'], d['T'], d['G']
    W = 2 * p.bit_length() + 4
    o = ['(set-logic QF_BV)']
    decls, asserts = [], []
    def C(x): return '(_ bv%d %d)' % (x % (1 << W), W)
    def V(name):
        decls.append('(declare-fun %s () (_ BitVec %d))' % (name, W))
        asserts.append('(assert (bvult %s %s))' % (name, C(p)))
        return name
    def mulmod(a, b, name):
        r = V(name); q = V(name + '_q')
        asserts.append('(assert (= (bvmul %s %s) (bvadd (bvmul %s %s) %s)))'
                       % (a, b, C(p), q, r))
        return r
    def red(expr, name):
        """expr is a BV term whose true value is in (-4p, 4p); reduce to [0,p)."""
        r = V(name)
        q = V(name + '_r')
        asserts.append('(assert (bvult %s %s))' % (q, C(8)))
        asserts.append('(assert (= (bvadd %s %s) (bvadd (bvmul %s %s) %s)))'
                       % (expr, C(4 * p), C(p), q, r))
        return r
    ax, ay = C(G[0]), C(G[1])
    for i in range(1, n):
        bx, by = lad[i]
        bi = 'b%d' % i
        decls.append('(declare-fun %s () (_ BitVec %d))' % (bi, W))
        asserts.append('(assert (bvult %s %s))' % (bi, C(2)))
        dd = red('(bvsub %s %s)' % (C(bx), ax), 'd%d' % i)
        nn = red('(bvsub %s %s)' % (C(by), ay), 'n%d' % i)
        asserts.append('(assert (not (= %s %s)))' % (dd, C(0)))
        d2 = mulmod(dd, dd, 'd2_%d' % i)
        n2 = mulmod(nn, nn, 'n2_%d' % i)
        sx = V('sx%d' % i); sy = V('sy%d' % i)
        t = red('(bvadd %s %s %s)' % (sx, ax, C(bx)), 't%d' % i)
        lhs = mulmod(t, d2, 'td%d' % i)
        asserts.append('(assert (= %s %s))' % (lhs, n2))
        u = red('(bvadd %s %s)' % (sy, ay), 'u%d' % i)
        v = red('(bvsub %s %s)' % (ax, sx), 'v%d' % i)
        l1 = mulmod(u, dd, 'ud%d' % i); l2 = mulmod(nn, v, 'nv%d' % i)
        asserts.append('(assert (= %s %s))' % (l1, l2))
        nax = V('ax%d' % i); nay = V('ay%d' % i)
        gx = mulmod(bi, red('(bvsub %s %s)' % (sx, ax), 'dx%d' % i), 'gx%d' % i)
        gy = mulmod(bi, red('(bvsub %s %s)' % (sy, ay), 'dy%d' % i), 'gy%d' % i)
        asserts.append('(assert (= %s %s))' % (nax, red('(bvadd %s %s)' % (ax, gx), 'mx%d' % i)))
        asserts.append('(assert (= %s %s))' % (nay, red('(bvadd %s %s)' % (ay, gy), 'my%d' % i)))
        ax, ay = nax, nay
    asserts.append('(assert (= %s %s))' % (ax, C(T[0])))
    asserts.append('(assert (= %s %s))' % (ay, C(T[1])))
    o += decls + asserts + ['(check-sat)']
    open(path, 'w').write('\n'.join(o) + '\n')
    return path
