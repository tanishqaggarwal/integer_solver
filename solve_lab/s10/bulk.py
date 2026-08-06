"""S11 step 42: BULK activation -- does the closure gain a kernel?

The closure is 1655 x 707 with rank 707: zero kernel, so the free inputs are
completely pinned.  Activation adds columns (s10/second.py: +1-2 knobs per pair).
Activate many dead paths at once and rebuild: if columns outgrow rank, a kernel
appears and the obstruction may become resolvable.
"""
import os, sys, random, collections, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE); FORBID = {2081, 4287}
random.seed(31)
v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
BAD = [21617, 29539]
CHECKS = sorted(a for a in range(L.NA) if a not in atom_out)

def cone(seeds):
    c, st = set(), list(seeds)
    while st:
        t = st.pop()
        if t in c: continue
        c.add(t)
        a = definer.get(t)
        if a is None: continue
        for w in L.avars[a]:
            if w != t: st.append(w)
    return c
seeds = set()
for a in BAD: seeds |= set(L.avars[a])
CC = cone(seeds)
pool = sorted(u for u in CC if u in FREE and v0[u] == 0 and u not in FORBID)
print(f'activation pool (zero free inputs in the cluster cone): {len(pool)}', flush=True)

def closure_of(v):
    vm = [x % P for x in v]
    gc = {}
    def gr(c):
        if c not in gc: gc[c] = set(ad.grad(c, vm)) - FORBID
        return gc[c]
    U, rows, cols = set(), set(BAD), {}
    for it in range(7):
        newU = set()
        for c in rows: newU |= gr(c)
        newU -= U; U |= newU
        for u in sorted(newU): cols[u] = jac_column(u, v, vm, CHECKS)
        nr = set(BAD)
        for u in U: nr |= set(cols[u])
        grew = nr - rows; rows = nr
        if not newU and not grew: break
    return sorted(rows), sorted(U), cols

def solve(rows, Us, cols, av):
    ri = {c: i for i, c in enumerate(rows)}
    n, m = len(rows), len(Us)
    M = [[0] * (m + 1) for _ in rows]
    for j, u in enumerate(Us):
        for c, d in cols[u].items():
            if c in ri: M[ri[c]][j] = d % P
    for c in rows: M[ri[c]][m] = (-av[c]) % P
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
    bad = sum(1 for i in range(r_, n) if M[i][m])
    return n, m, r_, bad

t0 = time.time()
for N in (0, 10, 30, 60, len(pool)):
    v = list(v0)
    for u in pool[:N]: v[u] = random.randrange(1, 1 << 48)
    ad.fwd(v, rounds=6)
    av = L.all_atom_values(v)
    nz = [a for a in range(L.NA) if av[a]]
    rows, Us, cols = closure_of(v)
    n, m, rank, bad = solve(rows, Us, cols, av)
    print(f'N={N:>4}: nonzero atoms {len(nz):>4}  closure {n} x {m}  rank {rank}  '
          f'kernel {m-rank}  inconsistent rows {bad}  score '
          f'{L.NEQ-len(L.failing_eqs(av))}  ({time.time()-t0:.0f}s)', flush=True)
