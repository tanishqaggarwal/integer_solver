#!/usr/bin/env python3
"""K31: the B half is derived FORWARD from its 78 selectors (K29 proved that), yet its value
is not the composition of the two live leaves.  So identify what the value actually IS."""
import sys, os, json, itertools
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from k26_drive import drive, rootpair, C, P

ch = json.load(open(K + '/chain.json'))
D = FD.points()
bypow = {}
for i_s, e in ch['exp'].items():
    bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
exp2sel = {ch['exp'][str(i)]: ch['sel'][str(i)] for i in range(256)}
S = FD.SHIFT


def neg(p): return (p[0], (-p[1]) % P)


def mul(m, Pt):
    R = FD.INF; Q = Pt; mm = abs(m)
    while mm:
        if mm & 1: R = FD.add(R, Q)
        Q = FD.add(Q, Q); mm >>= 1
    return R if m >= 0 else neg(R)


for pair in ([3, 5], [3, 10], [5, 10]):
    e1, e2 = pair
    v = drive({exp2sel[e1], exp2sel[e2]})
    A, B = rootpair(v)
    P1, P2 = bypow[e1], bypow[e2]
    print('\n=== ON = exponents %s ===' % pair)
    print('B =', B)
    print('  on cubic:', (B[1] ** 2 - pow(B[0], 3, P) - FD.B) % P == 0)
    # candidate 1: raw coordinate-wise sum (what a mux does if both leaves share an OR-group)
    cw = ((P1[0] + P2[0]) % P, (P1[1] + P2[1]) % P)
    cwS = ((P1[0] + P2[0] - S) % P, (P1[1] + P2[1]) % P)   # shifted-coord version
    print('  == coordwise sum (shifted coords):', B == cw, '  (raw x, minus one shift):', B == cwS)
    # candidate 2: composition, both signs
    print('  == P1+P2:', B == FD.add(P1, P2), '  == P1-P2:', B == FD.add(P1, neg(P2)),
          '  == P2-P1:', B == FD.add(P2, neg(P1)), '  == -(P1+P2):', B == neg(FD.add(P1, P2)))
    print('  == P1:', B == P1, '  == P2:', B == P2, '  == 2P1:', B == FD.add(P1, P1),
          '  == 2P2:', B == FD.add(P2, P2))
    # candidate 3: small multiples a*P1 + b*P2
    hit = None
    for a in range(-4, 5):
        for b in range(-4, 5):
            if a == 0 and b == 0: continue
            if FD.add(mul(a, P1), mul(b, P2)) == B: hit = (a, b)
    print('  small combo a*P1+b*P2 ==B :', hit)
    # candidate 4: is B some OTHER single leaf?
    lf = {bypow[e]: e for e in range(256)}
    print('  B is leaf exponent:', lf.get(B))
