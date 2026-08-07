#!/usr/bin/env python3
"""Plant a known low-weight answer as a FAKE complement target, so the negative result
can be trusted.  usage: yplant.py i,j,k,... out.txt"""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'ydata.json')))
p = int(d['p']); A_ = int(d['a']); N = int(d['N'])
L = [(int(a), int(b)) for a, b in d['ladder']]
def inv(z): return pow(z, p - 2, p)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if (x1 - x2) % p == 0:
        if (y1 + y2) % p == 0: return None
        l = (3 * x1 * x1 + A_) % p * inv(2 * y1 % p) % p
    else:
        l = (y2 - y1) % p * inv((x2 - x1) % p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)
def mul(kk, P):
    R = None; Q = P
    while kk > 0:
        if kk & 1: R = add(R, Q)
        Q = add(Q, Q); kk >>= 1
    return R
S = sorted(int(v) for v in sys.argv[1].split(','))
P = None
for i in S: P = add(P, L[i])
k = sum(1 << i for i in S)
assert mul(k % N, L[0]) == P
with open(sys.argv[2], 'w') as f:
    f.write('%d %d\n' % (P[0], P[1]))
    for x, y in L: f.write('%d %d\n' % (x, y))
print('planted S=%s |S|=%d  k=%d  -> %s' % (S, len(S), k, sys.argv[2]))
print('  (this S is the COMPLEMENT set; the corresponding forward weight is %d)' % (256 - len(S)))
