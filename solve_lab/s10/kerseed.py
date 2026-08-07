"""S11 step 45: the kernel of the compensation closure, and its support on the seed.

The closure is 433 x 454 -- more atoms than equations -- so a kernel exists.  If it
contains a vector nonzero on the seven residual atoms, then AT THE ATOM LEVEL every
one of those equations can hold while the seven are nonzero.  Compute it, find the
minimal support, and list exactly which extra atoms must go nonzero.
"""
import os, sys, collections, time, random
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
SEED = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
forced = set()
for e in range(L.NEQ):
    m, sq, co = L.eq_atoms[e]
    nz = [a for a, c in co.items() if c]
    if len(nz) == 1: forced.add(nz[0])

ACTIVE = set(SEED)
for rnd in range(14):
    OBLIG = set()
    for a in ACTIVE: OBLIG |= set(L.atom2eq[a])
    cand = set()
    for e in OBLIG:
        m, sq, co = L.eq_atoms[e]
        for a, c in co.items():
            if c and a not in forced: cand.add(a)
    new = cand - ACTIVE
    if not new: break
    ACTIVE |= new
OBLIG = set()
for a in ACTIVE: OBLIG |= set(L.atom2eq[a])
rows = sorted(OBLIG); cols = sorted(ACTIVE)
ci = {a: j for j, a in enumerate(cols)}
n, m = len(rows), len(cols)
print(f'closure {n} equations x {m} atoms', flush=True)

M = [[0] * m for _ in rows]
for i, e in enumerate(rows):
    mm, sq, co = L.eq_atoms[e]
    for a, c in co.items():
        if a in ci: M[i][ci[a]] = c % P
# eliminate mod p
piv, r_ = [], 0
t0 = time.time()
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
    piv.append(j); r_ += 1
print(f'rank {r_}; kernel dimension {m - r_}  ({time.time()-t0:.0f}s)', flush=True)
ps = set(piv)
free_cols = [j for j in range(m) if j not in ps]
def kervec(fc):
    z = [0] * m; z[fc] = 1
    for i, pj in enumerate(piv): z[pj] = (-M[i][fc]) % P
    return z
seedj = [ci[a] for a in SEED]
print(f'\nseed columns are pivots? {[cols[j] for j in seedj if j in ps]}')
print(f'seed columns that are FREE: {[cols[j] for j in seedj if j not in ps]}')
hits = []
for fc in free_cols:
    z = kervec(fc)
    s = [cols[j] for j in seedj if z[j]]
    if s: hits.append((len(s), cols[fc], s, sum(1 for x in z if x)))
hits.sort(reverse=True)
print(f'\nkernel basis vectors touching the seed: {len(hits)} of {len(free_cols)}')
for k, fc, s, supp in hits[:8]:
    print(f'  free col a{fc}: touches seed atoms {s}; total support {supp}')
if not hits:
    print('  NONE -- every kernel vector vanishes on all seven seed atoms')
    print('  => at the atom level the seven CANNOT all be compensated')
