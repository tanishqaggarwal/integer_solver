"""S11 step 66: separate the EXACTLY-LINEAR checks, then test that subsystem.

Every closure so far linearised at a point, and 656 of 1,376 large-move predictions
were wrong -- so those vetoes are untrustworthy.  But some checks ARE exactly linear
mod p in the free inputs (the cluster residues, verified).  Restricted to those, the
system is not a linearisation at all: it is exact.  If that subsystem is already
inconsistent, no assignment can satisfy it, and the instance is infeasible.

Usage: exactlin.py NTEST   (NTEST large random probes per free input)
"""
import os, sys, random, time, collections
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
random.seed(20260806)
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 3
v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
CHECKS = sorted(a for a in range(L.NA) if a not in L.atom_out)
BAD = [21617, 29539]
FORBID = {2081, 4287}
U = sorted((set(ad.grad(BAD[0], vm0)) | set(ad.grad(BAD[1], vm0))) - FORBID)
print(f'probing {len(U)} free inputs with {NT} large random moves each', flush=True)
cols = {u: jac_column(u, v0, vm0, CHECKS) for u in U}
lin_ok = collections.Counter()
lin_bad = collections.Counter()
t0 = time.time()
for i, u in enumerate(U):
    col = cols[u]
    for _ in range(NT):
        d = random.randrange(1, P)
        w = list(v0); w[u] = w[u] + d
        ad.fwd(w, rounds=6)
        aw = L.all_atom_values(w)
        for c in col:
            pred = (av0[c] + col[c] * d) % P
            if pred == aw[c] % P: lin_ok[c] += 1
            else: lin_bad[c] += 1
    if i % 20 == 0:
        print(f'  {i}/{len(U)}  ({time.time()-t0:.0f}s)', flush=True)
touched = set(lin_ok) | set(lin_bad)
exact = sorted(c for c in touched if lin_bad[c] == 0)
print(f'\nchecks touched: {len(touched)};  EXACTLY LINEAR in every probe: {len(exact)}')
print(f'  the two targets exactly linear? '
      f'{[c in exact for c in BAD]}')
# build the exact subsystem and test consistency
rows = exact
ri = {c: i for i, c in enumerate(rows)}
n, m = len(rows), len(U)
M = [[cols[u].get(c, 0) % P for u in U] for c in rows]
b = [(-av0[c]) % P for c in rows]
A = [M[i][:] + [b[i]] for i in range(n)]
piv, r_ = [], 0
for j in range(m):
    k = next((i for i in range(r_, n) if A[i][j]), None)
    if k is None: continue
    A[r_], A[k] = A[k], A[r_]
    inv = pow(A[r_][j], -1, P)
    A[r_] = [x * inv % P for x in A[r_]]
    for i in range(n):
        if i != r_ and A[i][j]:
            f = A[i][j]
            A[i] = [(x - f * z) % P for x, z in zip(A[i], A[r_])]
    piv.append(j); r_ += 1
bad = [i for i in range(r_, n) if A[i][m]]
print(f'\nEXACT subsystem: {n} rows x {m} cols, rank {r_}, inconsistent rows {len(bad)}')
if bad:
    print('*** the exactly-linear subsystem is ALREADY INCONSISTENT')
    print('    -> no assignment of these free inputs satisfies them all, and this')
    print('       is an exact statement, not a linearisation')
else:
    print('    the exactly-linear subsystem is CONSISTENT -- the obstruction lives')
    print('    entirely in the nonlinear checks, where no veto can be trusted')
