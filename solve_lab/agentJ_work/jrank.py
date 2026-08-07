#!/usr/bin/env python3
"""Bipartite structure of (equations x atoms) and rank of the coefficient matrix.

The parse says: eq_i  ==  mult_i * (sum_j c_ij A_j)^k_i  = 0, mult_i != 0.
Hence the whole system is equivalent to the LINEAR system  M a = 0  in the
39033 atom values a_j.  If rank(M) = #atoms then every atom must vanish.
"""
import pickle, os, sys, time, random
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
M = pickle.load(open(os.path.join(HERE, 'jmodel2.pkl'), 'rb'))
eqs, atoms = M['eqs'], M['atoms']
NA = len(atoms)

print("mult==0 count:", sum(1 for e in eqs if e['mult'] == 0))
colcnt = Counter()
for e in eqs:
    for c, j in e['terms']:
        colcnt[j] += 1
print("atoms not appearing:", NA - len(colcnt))
print("col count histogram:", Counter(colcnt.values()).most_common(10))

# connected components over the bipartite graph
parent = list(range(NA))
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: parent[ra] = rb

for e in eqs:
    ids = [j for _, j in e['terms']]
    for j in ids[1:]:
        union(ids[0], j)

comp = defaultdict(list)
for j in range(NA):
    comp[find(j)].append(j)
sizes = sorted((len(v) for v in comp.values()), reverse=True)
print("n components:", len(comp))
print("largest sizes:", sizes[:20])
print("size histogram:", Counter(sizes).most_common(10))

pickle.dump({'comp': {k: v for k, v in comp.items()}}, open(os.path.join(HERE, 'jcomp.pkl'), 'wb'))
