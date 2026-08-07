#!/usr/bin/env python3
"""Exhaustive low-Hamming-weight search for k with k*G = T (weight <= 6, unsigned)."""
import json, sys, time, itertools
from model import P, TARGET, to_short
from group import add, neg, mul
from fastgrp import batch_add

lad = json.load(open('ladder.json'))['ladder']
L = [(int(x), int(y)) for _, x, y in lad]
G = L[0]; T = to_short(TARGET)
assert mul(1 << 7, G) == L[7]
MASK = (1 << 62) - 1
t0 = time.time()

store = {}                     # trunc(x) -> k   for all weight<=3 sums
def put(pt, k):
    if pt is None: return
    store.setdefault(pt[0] & MASK, []).append(k)

for i in range(256): put(L[i], 1 << i)
pairs, pk = [], []
for i in range(256):
    for j in range(i + 1, 256):
        pairs.append((L[i], L[j])); pk.append((1 << i) | (1 << j))
res = batch_add([a for a, b in pairs], [b for a, b in pairs])
lvl2 = []
for r, k in zip(res, pk):
    if isinstance(r, tuple) and r and r[0] == 'EXC': continue
    put(r, k); lvl2.append((r, k))
print('weight<=2 stored, %d, %.1fs' % (len(store), time.time() - t0), flush=True)

CH = 200000
buf_p, buf_q, buf_k = [], [], []
n3 = 0
def flush():
    global buf_p, buf_q, buf_k, n3
    if not buf_p: return
    r = batch_add(buf_p, buf_q)
    for x, k in zip(r, buf_k):
        if isinstance(x, tuple) and x and x[0] == 'EXC': continue
        put(x, k); n3 += 1
    buf_p, buf_q, buf_k = [], [], []

for (pt, k) in lvl2:
    hi = k.bit_length()
    for m in range(hi, 256):
        buf_p.append(pt); buf_q.append(L[m]); buf_k.append(k | (1 << m))
        if len(buf_p) >= CH: flush()
flush()
print('weight<=3 stored: %d keys, %d level-3, %.1fs' % (len(store), n3, time.time() - t0), flush=True)

# query side: for every weight<=3 subset S, is T - S in the store?
hits = []
def check(pt, k):
    c = store.get(pt[0] & MASK)
    if not c: return
    for k2 in c:
        if mul(k + k2, G) == T:
            hits.append(k + k2); print('HIT k =', k + k2, flush=True)

nT = neg(T)
qb_p, qb_k = [], []
def qflush():
    global qb_p, qb_k
    if not qb_p: return
    r = batch_add([T] * len(qb_p), [neg(x) for x in qb_p])
    for x, k in zip(r, qb_k):
        if isinstance(x, tuple) and x and x[0] == 'EXC': continue
        if x is not None: check(x, k)
    qb_p, qb_k = [], []

# weight 0..3 on the query side
if T[0] & MASK in store: check(T, 0)
for i in range(256):
    qb_p.append(L[i]); qb_k.append(1 << i)
qflush()
cnt = 0
for (pt, k) in lvl2:
    qb_p.append(pt); qb_k.append(k)
    if len(qb_p) >= CH: qflush()
qflush()
print('query weight<=2 done %.1fs' % (time.time() - t0), flush=True)
for (pt, k) in lvl2:
    hi = k.bit_length()
    for m in range(hi, 256):
        r = None
        qb_p.append(pt); qb_k.append(k)   # placeholder replaced below
        qb_p.pop(); qb_k.pop()
# level-3 query: recompute level-3 points in chunks and query
buf_p, buf_q, buf_k = [], [], []
def flush_q():
    global buf_p, buf_q, buf_k
    if not buf_p: return
    r = batch_add(buf_p, buf_q)
    pts, ks = [], []
    for x, k in zip(r, buf_k):
        if isinstance(x, tuple) and x and x[0] == 'EXC': continue
        pts.append(x); ks.append(k)
    r2 = batch_add([T] * len(pts), [neg(x) for x in pts])
    for x, k in zip(r2, ks):
        if isinstance(x, tuple) and x and x[0] == 'EXC': continue
        if x is not None: check(x, k)
    buf_p, buf_q, buf_k = [], [], []
done = 0
for (pt, k) in lvl2:
    hi = k.bit_length()
    for m in range(hi, 256):
        buf_p.append(pt); buf_q.append(L[m]); buf_k.append(k | (1 << m))
        if len(buf_p) >= CH:
            flush_q(); done += CH
            print('  q %d  %.0fs' % (done, time.time() - t0), flush=True)
flush_q()
print('DONE weight<=6 exhaustive.  hits: %s   %.1fs' % (hits, time.time() - t0), flush=True)
json.dump({'hits': [str(h) for h in hits]}, open('runs/lowweight6.json', 'w'))
