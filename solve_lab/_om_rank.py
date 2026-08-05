#!/usr/bin/env python3
"""Sparse rank of the full 39033 x 39033 atom-coefficient matrix mod a prime,
using Markowitz (min row-count x col-count) pivoting."""
import pickle, sys, heapq, random
from collections import defaultdict

P = (1 << 31) - 1

D = pickle.load(open('_om_parsed2.pkl', 'rb'))
eqatoms = D['eqatoms']
atoms = sorted({k for d in eqatoms for k in d})
ai = {k: i for i, k in enumerate(atoms)}
N = len(atoms); M = len(eqatoms)
print('matrix %d x %d' % (M, N))

rows = [{ai[k]: c % P for k, c in d.items() if c % P} for d in eqatoms]
cols = defaultdict(set)
for r, d in enumerate(rows):
    for c in d: cols[c].add(r)

alive_row = [True] * M
rank = 0
# heap on (markowitz cost, col)
def colcost(c): return len(cols[c])
heap = [(len(cols[c]), c) for c in cols]
heapq.heapify(heap)
done_col = set()
t = 0
while heap:
    _, c = heapq.heappop(heap)
    if c in done_col: continue
    if not cols[c]:
        done_col.add(c); continue
    if _ != len(cols[c]):
        heapq.heappush(heap, (len(cols[c]), c)); continue
    # choose the row in cols[c] with fewest entries
    r = min(cols[c], key=lambda r: len(rows[r]))
    pr = rows[r]
    inv = pow(pr[c], P - 2, P)
    pr = {k: (v * inv) % P for k, v in pr.items()}
    rows[r] = pr
    targets = list(cols[c] - {r})
    for r2 in targets:
        f = rows[r2].get(c, 0)
        if not f: continue
        d2 = rows[r2]
        for k, v in pr.items():
            nv = (d2.get(k, 0) - f * v) % P
            if nv:
                if k not in d2: cols[k].add(r2)
                d2[k] = nv
            else:
                if k in d2:
                    del d2[k]; cols[k].discard(r2)
    for k in pr: cols[k].discard(r)
    alive_row[r] = False
    done_col.add(c)
    rank += 1
    t += 1
    if t % 2000 == 0:
        nnz = sum(len(rows[i]) for i in range(M) if alive_row[i])
        print('  eliminated %d  rank=%d  remaining nnz=%d' % (t, rank, nnz), flush=True)
    # refresh heap entries for touched cols
    for k in pr:
        if k not in done_col and cols[k]:
            heapq.heappush(heap, (len(cols[k]), k))

print('RANK =', rank, ' of N =', N, ' M =', M)
print('nullity (columns) =', N - rank)
leftover = [i for i in range(M) if alive_row[i] and rows[i]]
print('rows left with content:', len(leftover))
