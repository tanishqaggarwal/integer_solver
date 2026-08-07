#!/usr/bin/env python3
"""z3py builder for the ladder-fold feasibility problem (Int / BV), + SMT-LIB dump."""
from z3 import *

def build(d, logic='int', pin=None, W=None):
    p, n, lad, T, G = d['p'], d['n'], d['lad'], d['T'], d['G']
    C = []
    if logic == 'int':
        mk = lambda nm: Int(nm); K = lambda x: IntVal(x)
        lt = lambda a, b: a < b; ge = lambda a, b: a >= b
    else:
        W = W or 2 * p.bit_length() + 4
        mk = lambda nm: BitVec(nm, W); K = lambda x: BitVecVal(x % (1 << W), W)
        lt = lambda a, b: ULT(a, b); ge = lambda a, b: UGE(a, b)
    def V(nm):
        v = mk(nm); C.append(lt(v, K(p)))
        if logic == 'int': C.append(v >= 0)
        return v
    def mulmod(a, b, nm):
        r = V(nm); q = V(nm + '_q')
        C.append(a * b == K(p) * q + r)
        return r
    def red(e, nm, lo=4):
        r = V(nm); q = mk(nm + '_r')
        if logic == "int": C.append(And(q >= 0, q <= 2 * lo))
        else: C.append(ULT(q, K(2 * lo)))
        C.append(e + K(lo * p) == K(p) * q + r)
        return r
    ax, ay = K(G[0]), K(G[1])
    bits = []
    for i in range(1, n):
        bx, by = lad[i]
        bi = mk('b%d' % i); bits.append(bi)
        C.append(Or(bi == K(0), bi == K(1)))
        dd = red(K(bx) - ax, 'd%d' % i); nn = red(K(by) - ay, 'n%d' % i)
        C.append(dd != K(0))
        d2 = mulmod(dd, dd, 'd2_%d' % i); n2 = mulmod(nn, nn, 'n2_%d' % i)
        sx = V('sx%d' % i); sy = V('sy%d' % i)
        t = red(sx + ax + K(bx), 't%d' % i)
        C.append(mulmod(t, d2, 'td%d' % i) == n2)
        u = red(sy + ay, 'u%d' % i); vv = red(ax - sx, 'v%d' % i)
        C.append(mulmod(u, dd, 'ud%d' % i) == mulmod(nn, vv, 'nv%d' % i))
        nax = V('ax%d' % i); nay = V('ay%d' % i)
        gx = mulmod(bi, red(sx - ax, 'dx%d' % i), 'gx%d' % i)
        gy = mulmod(bi, red(sy - ay, 'dy%d' % i), 'gy%d' % i)
        C.append(nax == red(ax + gx, 'mx%d' % i))
        C.append(nay == red(ay + gy, 'my%d' % i))
        ax, ay = nax, nay
    C.append(ax == K(T[0])); C.append(ay == K(T[1]))
    if pin is not None:
        for i, bi in enumerate(bits, start=1):
            C.append(bi == K((pin >> i) & 1))
    return C, bits
