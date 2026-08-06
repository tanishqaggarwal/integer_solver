"""S11 step 22: cheapest sacrifice on the REDUCED (absorbable-veto) closure.

closure5.py steered the FULL closure and got 52 equations.  But the full closure
vetoes rows whose response to large moves is nonlinear, so its veto over-counts.
The reduced closure keeps only the p-absorbable rows, whose response is 90.6%
linear.  Force the cluster rows to be pivots, order the rest by decreasing cost,
and read off what must break.
"""
import os, sys, collections, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
P = ad.P
definer, atom_out = L.definer, L.atom_out
FREE = set(ad.FREE); FORBID = {2081, 4287}
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
CHECKS = set(a for a in range(L.NA) if a not in atom_out)
CH = sorted(CHECKS)
NOWFAIL = set(L.failing_eqs(av))
BAD = [21617, 29539]

def dz(a, w):
    s = 0
    for m, c in L.polys[a].items():
        k = m.count(w)
        if k == 0: continue
        if k == 1:
            t = c
            for z in m:
                if z != w: t *= v[z]
            s += t
        else: s += 2 * c * v[w]
    return s

absorb = collections.defaultdict(set)
for h in FREE:
    if h in FORBID: continue
    for a0 in L.var_atoms[h]:
        d = dz(a0, h)
        if d == 0 or d % P: continue
        if a0 in CHECKS: absorb[a0].add(h); continue
        t = atom_out[a0][1]
        if dz(a0, t) == 0: continue
        for c in L.var_atoms[t]:
            if c != a0 and c in CHECKS: absorb[c].add(h)
ABS = set(absorb)
gc = {}
def gr(c):
    if c not in gc: gc[c] = set(ad.grad(c, vm)) - FORBID
    return gc[c]
U, veto, cols = set(), set(BAD), {}
t0 = time.time()
for it in range(9):
    newU = set()
    for c in veto: newU |= gr(c)
    newU -= U; U |= newU
    for u in sorted(newU): cols[u] = jac_column(u, v, vm, CH)
    nv = set(BAD)
    for u in U: nv |= (set(cols[u]) & ABS)
    grew = nv - veto; veto = nv
    if not newU and not grew: break
print(f'reduced closure {len(veto)} x {len(U)}  ({time.time()-t0:.0f}s)', flush=True)

def cost1(c): return len(set(L.atom2eq[c]) - NOWFAIL)
Us = sorted(U); BADS = set(BAD)
order = sorted(veto, key=lambda c: (0 if c in BADS else 1, -cost1(c), c))
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
print(f'union of their equations: {len(eqs)}  -> score ceiling {L.NEQ - len(eqs)}')
print(f'(delivered 39026; canonical frame 39009)')
