#!/usr/bin/env python3
"""Backward cone of a set of variables through the syntactic definer map."""
import os, pickle, sys, json
from collections import defaultdict, deque, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
lead = pickle.load(open(os.path.join(HERE, 'jlead.pkl'), 'rb'))
atoms, polys = M['atoms'], P['polys']
NV = 38748
varsof = []
for p in polys:
    s = set()
    for k in p:
        s.update(k)
    varsof.append(s)
definer = {}
extra = defaultdict(list)
for i, v in enumerate(lead):
    if v is None:
        continue
    if v in definer:
        extra[v].append(i)
    else:
        definer[v] = i

d = json.load(open(sys.argv[1]))
val = [0] * NV
for k, v in d.items():
    val[int(k[2:]) if k.startswith('x_') else int(k)] = int(v)

start = [int(a) for a in sys.argv[2:]]
seen = set()
q = deque(start)
freev = []
levels = {}
for s in start:
    levels[s] = 0
while q:
    x = q.popleft()
    if x in seen:
        continue
    seen.add(x)
    i = definer.get(x)
    if i is None:
        freev.append(x)
        continue
    for w in varsof[i]:
        if w != x and w not in seen and w not in levels:
            levels[w] = levels[x] + 1
            q.append(w)
print("cone size:", len(seen), " free (no definer):", len(freev))
print("max depth:", max(levels.values()))
def show(x):
    s = str(x)
    return s if len(s) <= 26 else s[:14] + '..[%dd]' % len(s)
fs = sorted(freev, key=lambda z: levels[z])
print("\nfree inputs in cone (var, depth, value):")
for z in fs:
    print(f"  x_{z}  d={levels[z]}  {show(val[z])}")
