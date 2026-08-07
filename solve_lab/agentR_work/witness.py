#!/usr/bin/env python3
"""Python mirror of z3enc.build: computes the intended witness and checks every relation."""
def witness(d, k):
    p, n, lad, T, G = d['p'], d['n'], d['lad'], d['T'], d['G']
    W = {}
    def red(e, nm, lo=4):
        r = e % p; q = (e + lo * p - r) // p
        assert 0 <= r < p, (nm, r)
        assert 0 <= q <= 2 * lo, (nm, 'q', q, 'e', e)
        assert e + lo * p == p * q + r
        W[nm] = r; return r
    def mulmod(a, b, nm):
        r = a * b % p; q = (a * b - r) // p
        assert 0 <= r < p and 0 <= q < p, (nm, q)
        W[nm] = r; return r
    ax, ay = G
    for i in range(1, n):
        bx, by = lad[i]
        bi = (k >> i) & 1
        dd = red(bx - ax, 'd%d' % i); nn = red(by - ay, 'n%d' % i)
        assert dd != 0
        d2 = mulmod(dd, dd, 'd2_%d' % i); n2 = mulmod(nn, nn, 'n2_%d' % i)
        l = nn * pow(dd, p - 2, p) % p
        sx = (l * l - ax - bx) % p; sy = (l * (ax - sx) - ay) % p
        t = red(sx + ax + bx, 't%d' % i)
        assert mulmod(t, d2, 'td%d' % i) == n2, ('rel1', i)
        u = red(sy + ay, 'u%d' % i); vv = red(ax - sx, 'v%d' % i)
        assert mulmod(u, dd, 'ud%d' % i) == mulmod(nn, vv, 'nv%d' % i), ('rel2', i)
        gx = mulmod(bi, red(sx - ax, 'dx%d' % i), 'gx%d' % i)
        gy = mulmod(bi, red(sy - ay, 'dy%d' % i), 'gy%d' % i)
        nax = red(ax + gx, 'mx%d' % i); nay = red(ay + gy, 'my%d' % i)
        ax, ay = nax, nay
    assert (ax, ay) == T, ('target', (ax, ay), T)
    return W
