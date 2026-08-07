#!/usr/bin/env python3
"""Planted end-to-end test of the centre machinery at a NONTRIVIAL centre D.
Builds a fake T with a known S = D union E, then checks that the 'up' residual target
T - fold(D) is exactly fold(E) -- so the engine, unchanged, searches distance from D."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'ydata.json')))
p = int(d['p']); A_ = int(d['a'])
L = [(int(x), int(y)) for x, y in d['ladder']]
def inv(z): return pow(z, p - 2, p)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0: return None
        l = (3*x1*x1 + A_) % p * inv(2*y1 % p) % p
    else:
        l = (y2 - y1) % p * inv((x2 - x1) % p) % p
    x3 = (l*l - x1 - x2) % p
    return (x3, (l*(x1 - x3) - y1) % p)
def neg(P): return (P[0], (-P[1]) % p)
def fold(S):
    P = None
    for i in S: P = add(P, L[i])
    return P
D = [3, 17, 88, 201]
E = [11, 73, 140, 222, 244]
assert not set(D) & set(E)
S = sorted(D + E)
Tfake = fold(S)
up = add(Tfake, neg(fold(D)))
print('centre D            = %s' % D)
print('planted S = D u E   = %s   (|S| = %d, distance from D = %d)' % (S, len(S), len(E)))
print('up = T_fake - fold(D) == fold(E) : %s' % (up == fold(E)))
assert up == fold(E)
with open(os.path.join(HERE, 'plant_centre_up.txt'), 'w') as f:
    f.write('%d %d\n' % up)
    for x, y in L: f.write('%d %d\n' % (x, y))
print('wrote plant_centre_up.txt  (the engine should now find E, weight %d)' % len(E))
