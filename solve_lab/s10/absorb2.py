"""S11 step 7: proper p-absorbability, then the REDUCED system.

A check is p-absorbable if some free input h enters it with exact integer
derivative a nonzero multiple of p -- then h soaks up any multiple-of-p change and
the check needs only  residue == 0 (mod p).  Every other check must hold exactly,
and those are the nonlinear rows whose linear veto is untrustworthy.

Reduced system: rows = absorbable checks (mod-p preservation) + the two targets.
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
FREE = set(ad.FREE)
FORBID = {2081, 4287}
v = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm = [x % P for x in v]
av = L.all_atom_values(v)
CHECKS = set(a for a in range(L.NA) if a not in atom_out)

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
print(f'p-absorbable checks: {len(ABS)} of {len(CHECKS)}')
print(f'  a21617 {21617 in ABS} (handles {sorted(absorb.get(21617, []))[:4]})')
print(f'  a29539 {29539 in ABS} (handles {sorted(absorb.get(29539, []))[:4]})')

BAD = [21617, 29539]
U = sorted((set(ad.grad(BAD[0], vm)) | set(ad.grad(BAD[1], vm))) - FORBID)
cols = {u: jac_column(u, v, vm, sorted(CHECKS)) for u in U}
allrows = sorted(set().union(*[set(c) for c in cols.values()]))
rows = [c for c in allrows if c in ABS or c in BAD]
skipped = [c for c in allrows if c not in rows]
print(f'\naffected checks {len(allrows)}: {len(rows)} carry a mod-p veto, '
      f'{len(skipped)} are exact/nonlinear (veto dropped)')

ri = {c: i for i, c in enumerate(rows)}
n, m = len(rows), len(U)
M = [[0] * (m + 1) for _ in rows]
for j, u in enumerate(U):
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
print(f'reduced system {n} x {m}  rank {r_}  kernel {m-r_}  '
      f'inconsistent rows {len(bad_rows)}: {bad_rows[:8]}')
if not bad_rows:
    d = [0] * m
    for i, j in enumerate(piv): d[j] = M[i][m]
    sol = {U[j]: d[j] for j in range(m) if d[j]}
    print(f'*** CONSISTENT: delta moves {len(sol)} free inputs')
    json.dump({str(u): str(x) for u, x in sol.items()},
              open(os.path.join(HERE, 'delta_red.json'), 'w'))
    # also dump a kernel basis so we can search inside the solution space
    ker = []
    pivset = set(piv)
    for fc in range(m):
        if fc in pivset: continue
        z = [0] * m; z[fc] = 1
        for i, pj in enumerate(piv): z[pj] = (-M[i][fc]) % P
        ker.append(z)
    json.dump({'U': U, 'sol': [str(x) for x in d],
               'ker': [[str(x) for x in z] for z in ker]},
              open(os.path.join(HERE, 'redsol.json'), 'w'))
    print(f'    kernel dimension {len(ker)}; saved redsol.json')
