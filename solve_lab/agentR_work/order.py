#!/usr/bin/env python3
"""Group order of the fold law's group, by Cornacchia on 4p = L^2 + 27 M^2, verified on points."""
import math
from sympy import sqrt_mod, isprime, factorint
from model import P
from group import add, mul

def corn(p):
    for x0 in sqrt_mod(-27, p, all_roots=True):
        for x in (x0, p - x0, x0 + p, 2 * p - x0):
            if x % 2 == 0: continue
            a, b = 2 * p, x % (2 * p)
            lim = 2 * math.isqrt(p)
            while b * b > 4 * p:
                a, b = b, a % b
            c = 4 * p - b * b
            if c % 27: continue
            m2 = c // 27; m = math.isqrt(m2)
            if m * m == m2: return b, m
    return None

def orders(p):
    L, M = corn(p)
    ts = set()
    for s1 in (1, -1):
        ts.add(s1 * L)
        for s2 in (1, -1):
            for num in (s1 * L + 9 * s2 * M, s1 * L - 9 * s2 * M):
                if num % 2 == 0: ts.add(num // 2)
    return sorted({p + 1 - t for t in ts})

if __name__ == '__main__':
    import json
    from model import to_short, TARGET
    lad = json.load(open('ladder.json'))['ladder']
    G = (int(lad[0][1]), int(lad[0][2]))
    cands = orders(P)
    good = [n for n in cands if mul(n, G) is None]
    print('candidate orders:', len(cands))
    for n in good:
        print('order n =', n, ' prime?', isprime(n), ' bits', n.bit_length())
        print('   n-1 factors (partial):', factorint(n - 1, limit=10**6))
