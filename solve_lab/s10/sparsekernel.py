"""S10 step 65: is there a SPARSE vector in the wire kernel?

Cost of a kernel deformation = 0 identity equations + the square checks of the
members it moves.  A kernel vector moving only 2-3 members would cost ~3-6
equations and beat the 7 the current branch pays.

A kernel vector supported on S exists iff every kernel-coordinate row outside S
lies in a 2-dimensional subspace of Q^3.  So: view the kernel as a 220x3 matrix K
(row r = the three basis values at wire member r) and ask for the largest set of
rows lying in a common plane.  min|S| = (#nonzero rows) - (largest coplanar set).
"""
import os, sys, json, collections, itertools
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'wirekernel.json')))
WIRE, BASIS = d['wire'], d['basis']
N = len(WIRE)
K = [[BASIS[t][j] for t in range(3)] for j in range(N)]   # row per member

zero_rows = [j for j in range(N) if not any(K[j])]
nz_rows = [j for j in range(N) if any(K[j])]
print(f'wire members: {N}; kernel-zero rows: {len(zero_rows)} '
      f'({[WIRE[j] for j in zero_rows]}); nonzero rows: {len(nz_rows)}')

# group rows by direction (proportional rows) -- a rank-1 cluster
def norm(r):
    a, b, c = r
    g = 0
    for x in (a, b, c):
        g = __import__('math').gcd(g, abs(x))
    if g == 0:
        return (0, 0, 0)
    a, b, c = a // g, b // g, c // g
    for x in (a, b, c):
        if x != 0:
            if x < 0:
                a, b, c = -a, -b, -c
            break
    return (a, b, c)

dirs = collections.Counter(norm(K[j]) for j in nz_rows)
print(f'distinct directions among nonzero rows: {len(dirs)}')
print(f'largest rank-1 cluster: {dirs.most_common(1)[0][1]} rows')

# largest coplanar set: for each pair of independent directions, count rows in span
def det3(u, v, w):
    return (u[0]*(v[1]*w[2]-v[2]*w[1]) - u[1]*(v[0]*w[2]-v[2]*w[0])
            + u[2]*(v[0]*w[1]-v[1]*w[0]))

D = list(dirs)
print(f'searching planes over {len(D)} directions ...', flush=True)
best = (0, None)
for i in range(len(D)):
    for j in range(i + 1, len(D)):
        u, v = D[i], D[j]
        if det3(u, v, (1, 0, 0)) == 0 and det3(u, v, (0, 1, 0)) == 0 and det3(u, v, (0, 0, 1)) == 0:
            continue     # u,v parallel
        cnt = 0
        for dr, m in dirs.items():
            if det3(u, v, dr) == 0:
                cnt += m
        if cnt > best[0]:
            best = (cnt, (u, v))
print(f'largest coplanar set of nonzero rows: {best[0]} of {len(nz_rows)}')
minS = len(nz_rows) - best[0]
print(f'\n=> sparsest kernel vector has support >= {minS} members')
if minS <= 6:
    print('   *** SPARSE KERNEL VECTOR EXISTS -- wire route may come under budget')
    u, v = best[1]
    outside = [WIRE[j] for j in nz_rows if det3(u, v, norm(K[j])) != 0]
    print(f'   support: {outside}')
else:
    print('   no sparse kernel vector: every kernel deformation moves >= '
          f'{minS} members, so the wire route cannot get cheap this way')
