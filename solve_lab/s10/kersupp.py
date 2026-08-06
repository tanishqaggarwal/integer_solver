"""S11 step 46: what IS the support of a seed-touching kernel vector?

If those atoms can take the prescribed values, ALL their equations hold -- and the
set is closed, so that is a full solution.  Classify the support: how many are
p-absorbable gadgets (trivially settable through a free handle), how many are
gates, how many are pinned.
"""
import os, sys, collections, time
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
    OB = set()
    for a in ACTIVE: OB |= set(L.atom2eq[a])
    cand = set()
    for e in OB:
        m, sq, co = L.eq_atoms[e]
        for a, c in co.items():
            if c and a not in forced: cand.add(a)
    new = cand - ACTIVE
    if not new: break
    ACTIVE |= new
OB = set()
for a in ACTIVE: OB |= set(L.atom2eq[a])
rows = sorted(OB); cols = sorted(ACTIVE)
ci = {a: j for j, a in enumerate(cols)}
n, m = len(rows), len(cols)
M = [[0] * m for _ in rows]
for i, e in enumerate(rows):
    mm, sq, co = L.eq_atoms[e]
    for a, c in co.items():
        if a in ci: M[i][ci[a]] = c % P
piv, r_ = [], 0
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
ps = set(piv)
fc = ci[35758]                       # kernel vector touching a35758
z = [0] * m; z[fc] = 1
for i, pj in enumerate(piv): z[pj] = (-M[i][fc]) % P
supp = [cols[j] for j in range(m) if z[j]]
print(f'closure {n} x {m}; kernel vector on a35758 has support {len(supp)}')
print(f'  seed atoms in it: {[a for a in SEED if a in supp]}')

# classify the support
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
FREE = set(ad.FREE)
def dz(a, w):
    s = 0
    for mo, c in L.polys[a].items():
        k = mo.count(w)
        if k == 0: continue
        if k == 1:
            t = c
            for x in mo:
                if x != w: t *= v[x]
            s += t
        else: s += 2 * c * v[w]
    return s
kinds = collections.Counter()
absorbable = []
for a in supp:
    isgate = a in L.atom_out
    hasfree = any(w in FREE for w in set(L.avars[a]))
    # p-absorbable: a free input enters with derivative a nonzero multiple of p
    pabs = False
    for w in set(L.avars[a]):
        if w in FREE:
            d = dz(a, w)
            if d and d % P == 0: pabs = True; break
    kinds[('gate' if isgate else 'check', 'p-absorbable' if pabs else
           ('has-free-var' if hasfree else 'no-free-var'))] += 1
    if pabs: absorbable.append(a)
print(f'\nsupport composition: {dict(kinds)}')
print(f'  p-absorbable atoms in the support: {len(absorbable)} of {len(supp)}')
print(f'  currently nonzero among the support: '
      f'{[a for a in supp if av[a]]}')
eqs = set()
for a in supp: eqs |= set(L.atom2eq[a])
print(f'  the support touches {len(eqs)} equations (all inside the closure: '
      f'{eqs <= set(rows)})')
