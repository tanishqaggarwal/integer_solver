#!/usr/bin/env python3
"""Evaluate all atoms at a given assignment; report nonzero atoms and failing eqs."""
import os, pickle, sys, json
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
P = pickle.load(open(os.path.join(HERE, 'jpoly.pkl'), 'rb'))
eqs, atoms, polys = M['eqs'], M['atoms'], P['polys']
NV = 38748

path = sys.argv[1]
d = json.load(open(path))
v = [0] * NV
for k, val in d.items():
    v[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)

av = []
for p in polys:
    s = 0
    for k, c in p.items():
        t = c
        for i in k:
            t *= v[i]
        s += t
    av.append(s)
nz = [i for i, x in enumerate(av) if x]
print("nonzero atoms:", len(nz))
for i in nz[:40]:
    print(f"  atom {i}: {atoms[i][:90]}   value bits={av[i].bit_length()} val={av[i] if abs(av[i])<10**12 else str(av[i])[:40]+'...'}")

fails = []
for e in eqs:
    s = sum(c * av[j] for c, j in e['terms'])
    if s != 0:
        fails.append(e['i'])
print("failing eqs:", len(fails), fails[:20])
