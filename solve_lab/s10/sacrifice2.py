"""S10 step 107: exhaustive minimum sacrifice in the delivered frame (45 x 13)."""
import os, sys, json, itertools, collections, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import DETACH, definer, atom_out, ORDER, FREE, CHECKS, fwd, score, grad, jac_column
P = ad.P
w = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(w)
vm = [x % P for x in w]
av = L.all_atom_values(w)
BAD = [a for a in CHECKS if av[a]]
U = set()
for a in BAD: U |= set(grad(a, vm))
for c in (22230, 35762, 37887): U |= set(grad(c, vm))
U = sorted(U)
cols = {u: jac_column(u, vm) for u in U}
rows = sorted(set().union(*[set(cols[u]) for u in U]) | set(BAD))
print(f'frame-2 closure {len(rows)} rows x {len(U)} cols')
ri = {c: i for i, c in enumerate(rows)}

def consistent(drop):
    keep = [c for c in rows if c not in drop]
    n, m = len(keep), len(U)
    M = [[0] * (m + 1) for _ in keep]
    for i, c in enumerate(keep):
        for j, u in enumerate(U):
            d = cols[u].get(c)
            if d: M[i][j] = d % P
        M[i][m] = (-av[c]) % P
    r_ = 0
    for j in range(m):
        k = next((i for i in range(r_, n) if M[i][j]), None)
        if k is None: continue
        M[r_], M[k] = M[k], M[r_]
        inv = pow(M[r_][j], -1, P)
        M[r_] = [x * inv % P for x in M[r_]]
        for i in range(n):
            if i != r_ and M[i][j]:
                f = M[i][j]
                M[i] = [(a2 - f * b2) % P for a2, b2 in zip(M[i], M[r_])]
        r_ += 1
    return not any(M[i][m] for i in range(r_, n))

def cost(S):
    e = set()
    for c in S: e |= set(L.atom2eq[c])
    return len(e)

print(f'baseline (drop nothing) consistent? {consistent(set())}')
best = None
t0 = time.time()
for k in (1, 2, 3):
    found = []
    for S in itertools.combinations(rows, k):
        cst = cost(S)
        if best and cst >= best[0]: continue
        if consistent(set(S)):
            found.append((cst, S))
    found.sort()
    if found:
        print(f'size {k}: {len(found)} sets restore consistency; '
              f'cheapest {found[:5]}')
        if best is None or found[0][0] < best[0]: best = found[0]
    else:
        print(f'size {k}: none  ({time.time()-t0:.0f}s)')
if best:
    print(f'\n*** MINIMUM SACRIFICE in the delivered frame: {best[1]} '
          f'costing {best[0]} equations -> score {L.NEQ - best[0]}')
    for c in best[1]:
        print(f'   a{c}: {len(L.atom2eq[c])} eqs')
    json.dump({'drop': list(best[1]), 'cost': best[0]},
              open(os.path.join(HERE, 'sac2.json'), 'w'))
