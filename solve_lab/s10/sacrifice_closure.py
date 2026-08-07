"""S10 step 104: the exact minimum sacrifice for the TRUE closure (579 x 142).

Inconsistent  <=>  exists y in leftnull(J) with y.b != 0.
Dropping rows S restores consistency  <=>  t = Y.b lies in colspace(Y[:,S]).
Cost of dropping check c = number of equations containing atom c that are not
already failing.
"""
import os, sys, json, time, itertools
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
for a in BAD: U |= set(ad.grad(a, vm))
for c in [33796, 40562, 41400, 41507, 41827, 42245]: U |= set(ad.grad(c, vm))
U = sorted(U - FORBID)
cols = {u: jac_column(u, v, vm, CHECKS) for u in U}
rows = sorted(set().union(*[set(cols[u]) for u in U]) | set(BAD))
ri = {c: i for i, c in enumerate(rows)}
n, m = len(rows), len(U)
Jm = [[0] * m for _ in rows]
for j, u in enumerate(U):
    for c, d in cols[u].items(): Jm[ri[c]][j] = d % P
b = [(-av[c]) % P for c in rows]
print(f'closure {n} x {m}')

# leftnull(J) = null(J^T):  reduce J^T (m x n)
A = [[Jm[i][j] for i in range(n)] for j in range(m)]
piv, r_ = [], 0
for j in range(n):
    k = next((i for i in range(r_, m) if A[i][j]), None)
    if k is None: continue
    A[r_], A[k] = A[k], A[r_]
    inv = pow(A[r_][j], -1, P)
    A[r_] = [x * inv % P for x in A[r_]]
    for i in range(m):
        if i != r_ and A[i][j]:
            f = A[i][j]
            A[i] = [(x - f * y) % P for x, y in zip(A[i], A[r_])]
    piv.append(j); r_ += 1
rank = r_
free_cols = [j for j in range(n) if j not in set(piv)]
print(f'rank(J) = {rank};  leftnull dim = {len(free_cols)}')
Y = []                       # each row of Y is a left-null vector (length n)
for fc in free_cols:
    y = [0] * n
    y[fc] = 1
    for i, pj in enumerate(piv):
        y[pj] = (-A[i][fc]) % P
    Y.append(y)
K = len(Y)
t = [sum(Y[k][i] * b[i] for i in range(n)) % P for k in range(K)]
print(f't = Y.b  nonzero? {any(t)}  (nonzero entries {sum(1 for x in t if x)})')

def cost(c):
    eqs = set(L.atom2eq[c])
    now = set(L.failing_eqs(av))
    return len(eqs - now)

def in_span(S):
    basis = []
    for c in S:
        i = ri[c]
        w = [Y[k][i] for k in range(K)]
        for bp, bv in basis:
            if w[bp]:
                f = w[bp] * pow(bv[bp], -1, P) % P
                w = [(w[j] - f * bv[j]) % P for j in range(K)]
        nzj = next((j for j in range(K) if w[j]), None)
        if nzj is not None: basis.append((nzj, w))
    w = t[:]
    for bp, bv in basis:
        if w[bp]:
            f = w[bp] * pow(bv[bp], -1, P) % P
            w = [(w[j] - f * bv[j]) % P for j in range(K)]
    return not any(w)

cand = sorted(rows, key=cost)
print(f'\ncheapest rows: {[(c, cost(c)) for c in cand[:12]]}')
print('\n=== single-row sacrifices ===')
sols = [(cost(c), c) for c in rows if in_span([c])]
sols.sort()
if sols:
    for k, c in sols[:10]:
        print(f'  drop a{c}: RESTORES CONSISTENCY, cost {k} equations '
              f'-> score {L.NEQ - k}')
else:
    print('  none')
    print('\n=== pairs (cheapest 60 rows) ===')
    small = cand[:60]
    found = []
    for c1, c2 in itertools.combinations(small, 2):
        k = cost(c1) + cost(c2)
        if k >= 24: continue
        if in_span([c1, c2]): found.append((k, c1, c2))
    found.sort()
    print(f'  {len(found)} within-budget pairs restore consistency: {found[:8]}')
