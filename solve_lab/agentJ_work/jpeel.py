#!/usr/bin/env python3
"""Reverse-topological peel: repeatedly remove (var v, atom i) where v occurs in
only one surviving atom i and v is a def candidate of i.  Whatever survives is
the irreducible core."""
import os, pickle, time
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

occ = defaultdict(set)
for i, s in enumerate(varsof):
    for v in s:
        occ[v].add(i)

dcand = [set(d) for d in defcands]
alive = [True] * NA
order = []          # (atom, var) peel order, reverse-topological
queue = [v for v in occ if len(occ[v]) == 1]
t0 = time.time()
while queue:
    v = queue.pop()
    s = occ.get(v)
    if not s or len(s) != 1:
        continue
    i = next(iter(s))
    if not alive[i] or v not in dcand[i]:
        continue
    # peel
    alive[i] = False
    order.append((i, v))
    del occ[v]
    for w in varsof[i]:
        if w == v:
            continue
        o = occ.get(w)
        if o is None:
            continue
        o.discard(i)
        if len(o) == 1:
            queue.append(w)
        elif len(o) == 0:
            del occ[w]

print(f"peeled {len(order)} atoms/vars in {time.time()-t0:.1f}s")
print("surviving atoms:", sum(alive))
print("vars still occurring:", len(occ))
core_atoms = [i for i in range(NA) if alive[i]]
core_vars = sorted(occ)
print("core: atoms=%d vars=%d" % (len(core_atoms), len(core_vars)))
sg = Counter()
for i in core_atoms:
    p = polys[i]
    d = Counter(len(k) for k in p)
    sg[(d[0], d[1], d[2])] += 1
print("core atom signatures:", sg.most_common(20))
pickle.dump({'order': order, 'core_atoms': core_atoms, 'core_vars': core_vars},
            open(os.path.join(HERE, 'jpeel.pkl'), 'wb'))
