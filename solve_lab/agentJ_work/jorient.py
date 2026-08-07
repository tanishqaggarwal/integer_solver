#!/usr/bin/env python3
"""Efficient orientation: build a topological generative model.

Loop: fire every atom that has exactly one unknown variable which is a legal
definition target (linear, coefficient +-1, not in a higher monomial).  When
stuck, pick the atom with fewest unknowns and declare all-but-one of them FREE.
Result: a DAG order plus a set of free inputs plus a set of constraint atoms.
"""
import os, pickle, time, heapq, sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
polys, defcands = P['polys'], P['defcands']
NA = len(polys)
NV = 38748

varsof = []
for p in polys:
    s = set()
    for k in p:
        s.update(k)
    varsof.append(s)
occ = defaultdict(list)
for i, s in enumerate(varsof):
    for v in s:
        occ[v].append(i)
dset = [set(d) for d in defcands]

known = set()
unk = [len(varsof[i]) for i in range(NA)]
fired = [False] * NA
order = []
free = []
constraints = []
ready = [i for i in range(NA) if unk[i] == 1]
heap = []           # (unk_count, atom) lazy heap for stuck phase

def mark_known(v):
    known.add(v)
    for j in occ[v]:
        if not fired[j]:
            unk[j] -= 1
            if unk[j] == 1:
                ready.append(j)

def drain():
    while ready:
        i = ready.pop()
        if fired[i] or unk[i] != 1:
            continue
        miss = varsof[i] - known
        if len(miss) != 1:
            continue
        v = next(iter(miss))
        if v not in dset[i]:
            continue
        fired[i] = True
        order.append((i, v))
        mark_known(v)

t0 = time.time()
drain()
heap = [(unk[i], i) for i in range(NA) if not fired[i]]
heapq.heapify(heap)
while heap:
    u, i = heapq.heappop(heap)
    if fired[i]:
        continue
    miss = varsof[i] - known
    if len(miss) != unk[i] or unk[i] != u:
        if not fired[i]:
            unk[i] = len(miss)
            heapq.heappush(heap, (unk[i], i))
        continue
    if len(miss) == 0:
        fired[i] = True
        constraints.append(i)
        continue
    if len(miss) == 1 and next(iter(miss)) in dset[i]:
        v = next(iter(miss))
        fired[i] = True
        order.append((i, v)); mark_known(v); drain()
        continue
    # declare all but one free (prefer to keep a legal def target last)
    cand = [v for v in miss if v in dset[i]]
    keep = cand[0] if cand else None
    for v in miss:
        if v is keep:
            continue
        free.append(v); order.append((None, v)); mark_known(v)
    drain()
    if not fired[i]:
        miss = varsof[i] - known
        if len(miss) == 1 and next(iter(miss)) in dset[i]:
            v = next(iter(miss)); fired[i] = True
            order.append((i, v)); mark_known(v); drain()
        elif len(miss) == 0:
            fired[i] = True; constraints.append(i)
        else:
            for v in miss:
                free.append(v); order.append((None, v)); mark_known(v)
            fired[i] = True; constraints.append(i)
            drain()

print(f"done in {time.time()-t0:.1f}s")
print("free inputs:", len(free))
print("defined vars:", sum(1 for a, v in order if a is not None))
print("constraint atoms:", len(constraints))
print("known vars:", len(known), "of", NV)
pickle.dump({'order': order, 'free': free, 'constraints': constraints},
            open(os.path.join(HERE, 'jorient.pkl'), 'wb'))
