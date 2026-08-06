"""S10 step 112: the FULL closure -- expand columns from EVERY row, not just the
inconsistent ones.

closure2.py grew U only from the rows that witnessed the inconsistency.  But a
free input that touches any row of the system is a repair freedom, so the honest
fixed point is
    rows  <- checks moved by U
    U     <- free inputs that move any row
iterated together.  More columns can only help consistency.
"""
import os, sys, json, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
definer, atom_out = L.definer, L.atom_out
FORBID = {2081, 4287}
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
CHECKS = [a for a in range(L.NA) if a not in atom_out]
BAD = [a for a in CHECKS if av[a]]
print(f'failing {BAD}', flush=True)

gcache = {}
def gr(c):
    if c not in gcache: gcache[c] = set(ad.grad(c, vm)) - FORBID
    return gcache[c]

U, rows = set(), set(BAD)
cols = {}
t0 = time.time()
for it in range(8):
    newU = set()
    for c in rows: newU |= gr(c)
    newU -= U
    U |= newU
    for u in sorted(newU):
        cols[u] = jac_column(u, vm, CHECKS) if False else jac_column(u, v, vm, CHECKS)
    newrows = set()
    for u in U: newrows |= set(cols[u])
    newrows |= set(BAD)
    grew = newrows - rows
    rows = newrows
    print(f'it{it}: rows {len(rows)} (+{len(grew)})  cols {len(U)} (+{len(newU)})  '
          f'({time.time()-t0:.0f}s)', flush=True)
    if not newU and not grew: break

rows = sorted(rows); Us = sorted(U)
ri = {c: i for i, c in enumerate(rows)}
n, m = len(rows), len(Us)
print(f'\nFULL closure {n} x {m}; eliminating ...', flush=True)
M = [[0] * (m + 1) for _ in rows]
for j, u in enumerate(Us):
    for c, d in cols[u].items():
        if c in ri: M[ri[c]][j] = d % P
for c in rows: M[ri[c]][m] = (-av[c]) % P
rid = list(rows); piv, r_ = [], 0
for j in range(m):
    k = next((i for i in range(r_, n) if M[i][j]), None)
    if k is None: continue
    M[r_], M[k] = M[k], M[r_]; rid[r_], rid[k] = rid[k], rid[r_]
    inv = pow(M[r_][j], -1, P)
    M[r_] = [x * inv % P for x in M[r_]]
    for i in range(n):
        if i != r_ and M[i][j]:
            f = M[i][j]
            M[i] = [(a2 - f * b2) % P for a2, b2 in zip(M[i], M[r_])]
    piv.append(j); r_ += 1
bad_rows = [rid[i] for i in range(r_, n) if M[i][m]]
print(f'rank {r_} of {m} columns; inconsistent rows {len(bad_rows)}: {bad_rows[:10]}')
if not bad_rows:
    d = [0] * m
    for i, j in enumerate(piv): d[j] = M[i][m]
    sol = {Us[j]: d[j] for j in range(m) if d[j]}
    print(f'*** CONSISTENT: solution moves {len(sol)} free inputs')
    json.dump({str(u): str(x) for u, x in sol.items()},
              open(os.path.join(HERE, 'delta3.json'), 'w'))
    print('saved delta3.json')
else:
    print(f'kernel dim {m - r_}')
