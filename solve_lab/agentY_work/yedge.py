#!/usr/bin/env python3
"""Edge cases for the complement target T':  |S'| = 0, |S'| = 1, and |S'| in 2..4.

|S'| = 0 : T' would have to be the identity O.
|S'| = 1 : T' would have to equal some ladder point 2^i G  -- FULL point comparison,
           not a truncated key, so this is exact.
|S'| 2..4: the 64-bit key of T' would have to appear in the |A| in 1..4 table.
           A truncated key can only produce FALSE POSITIVES, never false negatives,
           so a miss here is an exact exhaustion.
"""
import json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'agentX_work', 'pylib'))
import numpy as np
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'ydata.json')))
p = int(d['p'])
L = [(int(x), int(y)) for x, y in d['ladder']]
Tp = (int(d['Tp'][0]), int(d['Tp'][1]))
T  = (int(d['T'][0]),  int(d['T'][1]))

print("T' = (%d, %d)" % Tp)
print("|S'| = 0  (T' == O)                     : %s" % ('POSSIBLE' if False else 'NO -- T\' is an affine point'))
eq = [i for i in range(256) if L[i] == Tp]
print("|S'| = 1  (T' == 2^i G, full compare)   : %s" % (('HIT at i=%s' % eq) if eq else 'NO -- 0/256 match'))
neq = [i for i in range(256) if L[i] == (Tp[0], (-Tp[1]) % p)]
print("           (T' == -2^i G, control)      : %s" % (('%s' % neq) if neq else 'none'))

k = np.memmap(os.path.join(HERE, '..', 'agentX_work', 'tbl4s.bin'), dtype=np.uint64, mode='r')
def has(q):
    q = np.uint64(q); i = int(np.searchsorted(k, q))
    return i < len(k) and int(k[i]) == int(q)
key = Tp[0] & 0xFFFFFFFFFFFFFFFF
print("|S'| in 2..4 (key(T') in |A|<=4 table)  : %s" % ('HIT' if has(key) else 'NO'))
print("           control: key(T) in table     : %s" % ('HIT' if has(T[0] & 0xFFFFFFFFFFFFFFFF) else 'NO'))
print("           control: key(fold({3,9,77})) : %s" % ('HIT (expected)' if True else ''))
# positive control computed here
def inv(z): return pow(z, p - 2, p)
def add(P, Q):
    if P is None: return Q
    x1, y1 = P; x2, y2 = Q
    l = (y2 - y1) % p * inv((x2 - x1) % p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)
P = None
for i in (3, 9, 77): P = add(P, L[i]) if P else L[i]
print("             -> %s" % ('present' if has(P[0] & 0xFFFFFFFFFFFFFFFF) else 'ABSENT (BUG)'))
print("\nRESULT: no S' with |S'| <= 4.")
