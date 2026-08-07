"""Classify the 10,792 check atoms by shape."""
import collections, sys
import dlib as L
P = L.P

# which free inputs occur in exactly one atom (handles)
occ = collections.Counter()
for a in range(L.NA):
    for u in L.avars[a]:
        occ[u] += 1
solo = {u for u in L.freeset if occ[u] == 1}
print('free inputs:', len(L.freeset), ' solo-occurrence free inputs (handles):', len(solo))

cnt = collections.Counter()
eqpins = []
for a in range(L.NA):
    if a in L.atom_out:
        continue
    p = L.polys[a]
    deg = max((len(m) for m in p), default=0)
    nt = len(p)
    cnt[(deg, nt)] += 1
print('check shape (deg, nterms) histogram, top 25:')
for k, n in cnt.most_common(25):
    print('  ', k, n)
