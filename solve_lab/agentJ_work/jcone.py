#!/usr/bin/env python3
"""Explore the neighbourhood of given variables: atoms they appear in, values."""
import os, pickle, sys, json
from collections import defaultdict, deque

HERE = os.path.dirname(os.path.abspath(__file__))
M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
atoms, polys = M['atoms'], P['polys']
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

d = json.load(open(sys.argv[1]))
v = [0] * NV
for k, val in d.items():
    v[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)

def av(i):
    s = 0
    for k, c in polys[i].items():
        t = c
        for j in k:
            t *= v[j]
        s += t
    return s

def show(x):
    if x == 0: return '0'
    s = str(x)
    return s if len(s) <= 24 else s[:12] + '..[%dd]' % len(s)

start = [int(a) for a in sys.argv[2:]]
seen = set()
q = deque(start)
depth = {x: 0 for x in start}
MAXD = int(os.environ.get('MAXD', '2'))
while q:
    x = q.popleft()
    if x in seen: continue
    seen.add(x)
    print(f"\n=== x_{x} = {show(v[x])}   (depth {depth[x]}, {len(occ[x])} atoms)")
    for i in occ[x]:
        print(f"    a{i}: {atoms[i][:120]}   -> {show(av(i))}")
        if depth[x] < MAXD:
            for w in varsof[i]:
                if w not in seen and w not in depth:
                    depth[w] = depth[x] + 1
                    q.append(w)
