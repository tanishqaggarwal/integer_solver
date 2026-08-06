"""S11 step 8: closure over ABSORBABLE rows only.

Veto rows = checks that can absorb a multiple of p through a free handle; those
genuinely must be preserved mod p.  Non-absorbable checks are dropped from the
veto: they are largely SHADOWS (exact integer multiples of a gadget, e.g.
a37662 = 10*a21617) which vanish automatically once their gadget does, and their
response to large moves is nonlinear anyway.  Expand columns from every veto row.
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
BAD = [21617, 29539]
gcache = {}
def gr(c):
    if c not in gcache: gcache[c] = set(ad.grad(c, vm)) - FORBID
    return gcache[c]

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
    print(f'it{it}: veto rows {len(veto)} (+{len(grew)})  cols {len(U)} (+{len(newU)})'
          f'  ({time.time()-t0:.0f}s)', flush=True)
    if not newU and not grew: break

rows = sorted(veto); Us = sorted(U)
ri = {c: i for i, c in enumerate(rows)}
n, m = len(rows), len(Us)
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
print(f'\nreduced closure {n} x {m}  rank {r_}  kernel {m-r_}  '
      f'inconsistent {len(bad_rows)}: {bad_rows[:10]}')
if not bad_rows:
    d = [0] * m
    for i, j in enumerate(piv): d[j] = M[i][m]
    pivset = set(piv)
    ker = []
    for fc in range(m):
        if fc in pivset: continue
        z = [0] * m; z[fc] = 1
        for i, pj in enumerate(piv): z[pj] = (-M[i][fc]) % P
        ker.append(z)
    json.dump({'U': Us, 'sol': [str(x) for x in d],
               'ker': [[str(x) for x in z] for z in ker]},
              open(os.path.join(HERE, 'redsol.json'), 'w'))
    print(f'*** CONSISTENT: delta moves {sum(1 for x in d if x)} inputs; '
          f'kernel dim {len(ker)}; saved redsol.json')
