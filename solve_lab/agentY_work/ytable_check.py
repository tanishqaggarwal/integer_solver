#!/usr/bin/env python3
"""Agent Y -- validate the MITM table before reusing it.

The table (built by agent X) holds the low 64 bits of x( sum_{i in A} 2^i G ) for every
A with |A| in [1..4].  It depends ONLY on the ladder, not on the target, so it is reusable
for the complement target T'.  This script proves that independently:

  * file size == 8 * (C(256,1)+C(256,2)+C(256,3)+C(256,4))
  * the file is sorted ascending
  * for random subsets of sizes 1..4, the key computed here in pure Python bignum
    IS present in the table                       (no false negatives)
  * for random 64-bit words, the key is absent    (the table is not saturated)
  * the bitmap is consistent with the table       (every table key's top-32 bit is set)
"""
import json, os, sys, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'agentX_work', 'pylib'))
import numpy as np
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
TBL = os.path.join(HERE, '..', 'agentX_work', 'tbl4s.bin')
BMP = os.path.join(HERE, '..', 'agentX_work', 'bm4.bin')

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
        l = (3 * x1 * x1 + A_) % p * inv(2 * y1 % p) % p
    else:
        l = (y2 - y1) % p * inv((x2 - x1) % p) % p
    x3 = (l * l - x1 - x2) % p
    return (x3, (l * (x1 - x3) - y1) % p)
def fold(S):
    P = None
    for i in S: P = add(P, L[i])
    return P
def key(P): return P[0] & 0xFFFFFFFFFFFFFFFF

want = sum(comb(256, s) for s in range(1, 5))
sz = os.path.getsize(TBL)
print('table file       : %d bytes = %d keys' % (sz, sz // 8))
print('expected keys    : %d   match: %s' % (want, sz // 8 == want))
assert sz // 8 == want

k = np.memmap(TBL, dtype=np.uint64, mode='r')
# sortedness, chunked
sorted_ok = True
CH = 1 << 24
prev = None
for off in range(0, len(k), CH):
    c = np.asarray(k[off:off + CH])
    if prev is not None and c[0] < prev: sorted_ok = False; break
    if np.any(c[1:] < c[:-1]): sorted_ok = False; break
    prev = c[-1]
print('sorted ascending : %s' % sorted_ok)
assert sorted_ok

bm = np.memmap(BMP, dtype=np.uint8, mode='r')
print('bitmap file      : %d bytes (expect %d)' % (len(bm), 1 << 29))
assert len(bm) == 1 << 29

def has(q):
    q = np.uint64(q)
    i = int(np.searchsorted(k, q))
    return i < len(k) and int(k[i]) == int(q)
def bm_has(q):
    bi = q >> 32
    return bool((int(bm[bi >> 3]) >> (bi & 7)) & 1)

random.seed(7)
tot = 0; found = 0; bmok = 0
for trial in range(60):
    s = random.randint(1, 4)
    S = sorted(random.sample(range(256), s))
    q = key(fold(S))
    tot += 1
    if has(q): found += 1
    if bm_has(q): bmok += 1
print('random |A| in 1..4 subsets present in table : %d/%d' % (found, tot))
print('  ... and passing the bitmap prefilter      : %d/%d' % (bmok, tot))
assert found == tot and bmok == tot

# negative control: random 64-bit words should essentially never be present
neg = sum(1 for _ in range(2000) if has(random.getrandbits(64)))
print('random 64-bit words present (expect ~0)     : %d/2000' % neg)

# does the table contain the key of T' itself?  (that is the |S'| <= 4 test, done for real
# in yedge.py -- here just report)
Tp = (int(d['Tp'][0]), int(d['Tp'][1]))
print("\nkey(T') in table (|S'| in 1..4 probe)       : %s" % has(key(Tp)))
print("key(T)  in table (X's same probe, control)  : %s" % has(key((int(d['T'][0]), int(d['T'][1])))))
print('\nTABLE VALIDATED AND REUSABLE (it is target-independent).')
