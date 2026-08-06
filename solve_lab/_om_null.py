#!/usr/bin/env python3
"""Null-space of the atom-coefficient matrix C[E(A),A] for a growing atom set A."""
import pickle, sys
from fractions import Fraction
from collections import defaultdict
import heal_harness as H
import _om_parse as OP

D = pickle.load(open('_om_parsed2.pkl', 'rb'))
eqatoms = D['eqatoms']; astof = D['astof']
ainc = defaultdict(list)
for e, d in enumerate(eqatoms):
    for k in d: ainc[k].append(e)
F = [2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]

def nullspace(M, nc):
    """exact rational nullspace basis of M (list of rows len nc); returns basis vectors"""
    M = [[Fraction(x) for x in r] for r in M]
    nr = len(M); piv = []
    r = 0
    for c in range(nc):
        s = None
        for i in range(r, nr):
            if M[i][c] != 0: s = i; break
        if s is None: continue
        M[r], M[s] = M[s], M[r]
        pv = M[r][c]
        M[r] = [x / pv for x in M[r]]
        for i in range(nr):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == nr: break
    free = [c for c in range(nc) if c not in piv]
    basis = []
    for fc in free:
        v = [Fraction(0)] * nc; v[fc] = Fraction(1)
        for i, c in enumerate(piv): v[c] = -M[i][fc]
        # clear denominators
        from math import gcd
        den = 1
        for x in v: den = den * x.denominator // gcd(den, x.denominator)
        vi = [int(x * den) for x in v]
        g = 0
        for x in vi: g = gcd(g, abs(x))
        if g > 1: vi = [x // g for x in vi]
        basis.append(vi)
    return basis, len(piv)

def analyze(A, label=''):
    A = sorted(A)
    E = set(F)
    for k in A: E |= set(ainc[k])
    E = sorted(E)
    M = [[eqatoms[e].get(k, 0) for k in A] for e in E]
    b, rk = nullspace(M, len(A))
    print('%s |A|=%d |E|=%d rank=%d nullity=%d' % (label, len(A), len(E), rk, len(b)))
    return A, E, b

if __name__ == '__main__':
    A0 = set()
    for e in F: A0 |= set(eqatoms[e])
    A, E, b = analyze(A0, 'F-atoms:')
    for v in b:
        nz = [(A[i], v[i]) for i in range(len(A)) if v[i] != 0]
        print('  nullvec (%d nz):' % len(nz), nz[:8])
