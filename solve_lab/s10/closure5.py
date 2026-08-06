"""S10 step 115: cheapest sacrifice on the FULL closure (1655 x 707).

Which rows end up witnessing the inconsistency depends entirely on pivot order.
Eliminate with rows sorted by DECREASING cost: expensive rows become pivots, so
whatever is left over -- and therefore whatever must be sacrificed -- is as cheap
as the greedy can make it.  Cost of a set = |union of its equations|.
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
NOWFAIL = set(L.failing_eqs(av))

gcache = {}
def gr(c):
    if c not in gcache: gcache[c] = set(ad.grad(c, vm)) - FORBID
    return gcache[c]

U, rows, cols = set(), set(BAD), {}
t0 = time.time()
for it in range(8):
    newU = set()
    for c in rows: newU |= gr(c)
    newU -= U; U |= newU
    for u in sorted(newU): cols[u] = jac_column(u, v, vm, CHECKS)
    nr = set(BAD)
    for u in U: nr |= set(cols[u])
    grew = nr - rows; rows = nr
    if not newU and not grew: break
print(f'closure {len(rows)} x {len(U)}  ({time.time()-t0:.0f}s)', flush=True)

def cost1(c): return len(set(L.atom2eq[c]) - NOWFAIL)
Us = sorted(U)
# expensive rows first  ->  cheap rows fall through as witnesses
BADS = set(BAD)
# force the failing rows to be PIVOTS (rank 0 in the order), then most expensive
order = sorted(rows, key=lambda c: (0 if c in BADS else 1, -cost1(c), c))
ri = {c: i for i, c in enumerate(order)}
n, m = len(order), len(Us)
M = [[0] * (m + 1) for _ in order]
for j, u in enumerate(Us):
    for c, d in cols[u].items():
        if c in ri: M[ri[c]][j] = d % P
for c in order: M[ri[c]][m] = (-av[c]) % P
rid = list(order); r_ = 0
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
    r_ += 1
wit = [rid[i] for i in range(r_, n) if M[i][m]]
eqs = set()
for c in wit: eqs |= set(L.atom2eq[c])
print(f'rank {r_}; witnesses ({len(wit)}): {[(c, len(L.atom2eq[c])) for c in wit]}')
print(f'union of their equations: {len(eqs)}  '
      f'(currently failing {len(NOWFAIL)})')
print(f'=> insisting on the cluster costs at least {len(eqs)} equations '
      f'-> score at most {L.NEQ - len(eqs)}')
print(f'   (current deliverable 39026; canonical frame now 39009)')
json.dump({'wit': wit, 'cost': len(eqs)}, open(os.path.join(HERE, 'wit5.json'), 'w'))
