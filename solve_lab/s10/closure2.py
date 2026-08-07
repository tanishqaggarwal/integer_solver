"""S10 step 103: closure loop on the linearised repair.

Solve J.delta == -r mod p.  If inconsistent, the elimination hands back rows that
witness the inconsistency; pull in every free input that can move THOSE rows and
re-solve.  Iterate to a fixed point.  This is the honest closure -- earlier
sessions closed only over the columns reachable from the failing checks.
"""
import os, sys, json, time
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

U = set()
for a in BAD:
    U |= set(ad.grad(a, vm))
U -= FORBID
cols = {}

def solve(U):
    Us = sorted(U)
    for u in Us:
        if u not in cols:
            cols[u] = jac_column(u, v, vm, CHECKS)
    rows = sorted(set().union(*[set(cols[u]) for u in Us]) | set(BAD))
    ri = {c: i for i, c in enumerate(rows)}
    n, m = len(rows), len(Us)
    M = [[0] * (m + 1) for _ in rows]
    for j, u in enumerate(Us):
        for c, d in cols[u].items():
            M[ri[c]][j] = d % P
    for c in rows:
        M[ri[c]][m] = (-av[c]) % P
    rid = list(rows)
    piv, r_ = [], 0
    for j in range(m):
        k = next((i for i in range(r_, n) if M[i][j]), None)
        if k is None: continue
        M[r_], M[k] = M[k], M[r_]; rid[r_], rid[k] = rid[k], rid[r_]
        inv = pow(M[r_][j], -1, P)
        M[r_] = [x * inv % P for x in M[r_]]
        for i in range(n):
            if i != r_ and M[i][j]:
                f = M[i][j]
                M[i] = [(a - f * b) % P for a, b in zip(M[i], M[r_])]
        piv.append(j); r_ += 1
    bad_rows = [rid[i] for i in range(r_, n) if M[i][m]]
    sol = None
    if not bad_rows:
        d = [0] * m
        for i, j in enumerate(piv): d[j] = M[i][m]
        sol = {Us[j]: d[j] for j in range(m) if d[j]}
    return rows, Us, r_, bad_rows, sol

t0 = time.time()
for it in range(12):
    rows, Us, rank, bad_rows, sol = solve(U)
    print(f'it{it}: {len(rows)} rows x {len(Us)} cols  rank {rank}  '
          f'inconsistent rows {len(bad_rows)}  ({time.time()-t0:.0f}s)', flush=True)
    if sol is not None:
        print(f'  *** CONSISTENT: solution moves {len(sol)} free inputs')
        json.dump({str(u): str(d) for u, d in sol.items()},
                  open(os.path.join(HERE, 'delta.json'), 'w'))
        print('  saved delta.json'); break
    new = set()
    for c in bad_rows[:40]:
        new |= set(ad.grad(c, vm))
    new -= FORBID | U
    print(f'  witnesses {bad_rows[:8]}; adding {len(new)} new free inputs', flush=True)
    if not new:
        print('  fixed point reached -- genuinely inconsistent over ALL reachable '
              'free inputs'); break
    U |= new
