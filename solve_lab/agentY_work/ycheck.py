#!/usr/bin/env python3
"""Cross-check the C engine's field/point arithmetic against Python bignum, for the
COMPLEMENT base T' (an independent check from X's, whose base was T)."""
import json, os, subprocess, sys
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'ydata.json')))
p = int(d['p']); A_ = int(d['a'])
L = [(int(x), int(y)) for x, y in d['ladder']]
Tp = (int(d['Tp'][0]), int(d['Tp'][1]))

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
def neg(P): return (P[0], (-P[1]) % p)

def limbs(v): return [(v >> (64 * i)) & 0xFFFFFFFFFFFFFFFF for i in range(4)]

out = subprocess.run([os.path.join(HERE, 'ymitm'), 'selftest',
                      os.path.join(HERE, 'data_comp.txt')],
                     capture_output=True, text=True).stdout
got = {}
for line in out.splitlines():
    if '=' in line and ('L5+L7' in line or 'T-L3' in line):
        lhs, rhs = line.split('=')
        got[lhs.strip()] = [int(t) for t in rhs.split()]

exp = {'L5+L7 x': limbs(add(L[5], L[7])[0]),
       'T-L3 x':  limbs(add(Tp, neg(L[3]))[0]),
       "T-L3 y":  limbs(add(Tp, neg(L[3]))[1])}
ok = True
for kname, v in exp.items():
    m = got.get(kname) == v
    ok &= m
    print('%-10s C==Python limb-for-limb : %s' % (kname, m))
print('\nARITHMETIC CROSS-CHECK: %s' % ('PASS' if ok else 'FAIL'))
sys.exit(0 if ok else 1)
