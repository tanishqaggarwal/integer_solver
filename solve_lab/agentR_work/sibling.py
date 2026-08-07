#!/usr/bin/env python3
"""Scaled-down siblings of the reduced problem: same law (y^2 = x^3 + b, chord addition),
same shape (ladder of 2^i G, boolean selectors, one target), but an m-bit prime instead of 256.
Used to measure how automated-reasoning solvers scale on THIS problem shape."""
import math, random
from sympy import isprime, sqrt_mod

def corn(p):
    for x0 in sqrt_mod(-27, p, all_roots=True):
        for x in (x0, p - x0):
            for xx in (x, x + p):
                if xx % 2 == 0: continue
                a, b = 2 * p, xx % (2 * p)
                while b * b > 4 * p: a, b = b, a % b
                c = 4 * p - b * b
                if c % 27: continue
                m2 = c // 27; mm = math.isqrt(m2)
                if mm * mm == m2: return b, mm
    return None

class Grp:
    def __init__(self, p, b):
        self.p, self.b = p, b
    def add(self, A, B):
        p = self.p
        if A is None: return B
        if B is None: return A
        ax, ay = A; bx, by = B
        if ax == bx:
            if (ay + by) % p == 0: return None
            l = 3 * ax * ax % p * pow(2 * ay, p - 2, p) % p
        else:
            l = (by - ay) % p * pow(bx - ax, p - 2, p) % p
        cx = (l * l - ax - bx) % p
        return (cx, (l * (ax - cx) - ay) % p)
    def mul(self, k, A):
        R = None; Q = A
        while k:
            if k & 1: R = self.add(R, Q)
            Q = self.add(Q, Q); k >>= 1
        return R

def make(mbits, seed=0):
    """prime p of mbits bits, p=1 mod 3, curve y^2=x^3+b of PRIME order N, generator G."""
    rnd = random.Random(1000 * mbits + seed)
    while True:
        p = rnd.randrange(1 << (mbits - 1), 1 << mbits) | 1
        if p % 3 != 1 or not isprime(p): continue
        c = corn(p)
        if c is None: continue
        L, M = c
        ts = {L, -L}
        for s in (1, -1):
            for u in (L + 9 * s * M, L - 9 * s * M, -L + 9 * s * M, -L - 9 * s * M):
                if u % 2 == 0: ts.add(u // 2)
        orders = sorted({p + 1 - t for t in ts})
        for b in range(1, 60):
            g = Grp(p, b)
            # a point on the curve
            G = None
            for x in range(1, min(p, 4000)):
                r = (x * x * x + b) % p
                y = sqrt_mod(r, p)
                if y is not None: G = (x, int(y)); break
            if G is None: continue
            ok = [N for N in orders if isprime(N) and g.mul(N, G) is None]
            if not ok: continue
            N = ok[0]
            return dict(p=p, b=b, N=N, G=G, grp=g)

def instance(mbits, nbits=None, seed=0):
    d = make(mbits, seed)
    g, N, G = d['grp'], d['N'], d['G']
    n = nbits or N.bit_length()
    rnd = random.Random(7 * mbits + 13 * seed + 3)
    k = rnd.randrange(1, 1 << n) | 1          # odd, so bit 0 is always on
    T = g.mul(k, G)
    lad = [G]
    for i in range(1, n): lad.append(g.add(lad[-1], lad[-1]))
    d.update(n=n, k=k, T=T, lad=lad)
    return d

if __name__ == '__main__':
    for m in (8, 12, 16, 20, 24, 28, 32):
        d = instance(m)
        print('m=%2d  p=%-12d b=%-3d N=%-12d n=%2d  k=%d' % (m, d['p'], d['b'], d['N'], d['n'], d['k']))
