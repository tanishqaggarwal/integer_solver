#!/usr/bin/env python3
"""Forward topological closure: repeatedly fire an atom that has exactly one
unknown variable which is a def candidate of it.  Reports the fired DAG, the
free inputs, and the leftover constraint atoms."""
import os, pickle, time, sys
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

def closure(seed):
    """seed: set of variables declared known.  Returns (order, known, fired)."""
    known = set(seed)
    unk = [len(varsof[i] - known) for i in range(NA)]
    fired = [False] * NA
    order = []
    q = [i for i in range(NA) if unk[i] == 1]
    while q:
        i = q.pop()
        if fired[i] or unk[i] != 1:
            continue
        miss = varsof[i] - known
        if len(miss) != 1:
            continue
        v = next(iter(miss))
        if v not in defcands[i]:
            continue
        fired[i] = True
        order.append((i, v))
        known.add(v)
        for j in occ[v]:
            if not fired[j]:
                unk[j] -= 1
                if unk[j] == 1:
                    q.append(j)
    return order, known, fired

t0 = time.time()
order, known, fired = closure(set())
print(f"closure from empty seed: fired {len(order)} atoms, known {len(known)} vars ({time.time()-t0:.1f}s)")
left = [i for i in range(NA) if not fired[i]]
print("unfired atoms:", len(left))
unknownvars = set(range(NV)) - known
print("unknown vars:", len(unknownvars))
sg = Counter()
for i in left:
    p = polys[i]
    d = Counter(len(k) for k in p)
    sg[(d[0], d[1], d[2])] += 1
print("unfired atom signatures:", sg.most_common(20))
# how many unknowns per unfired atom
print("unknowns/unfired atom:", Counter(len(varsof[i] - known) for i in left).most_common(12))
pickle.dump({'order': order, 'known': known, 'fired': fired},
            open(os.path.join(HERE, 'jfwd0.pkl'), 'wb'))
