"""S11 step 48: sparsify the settable-only kernel and classify its support.

An 8-dimensional kernel touching the seed with support ~71.  Greedily cancel
coordinates using the other basis vectors to find the sparsest seed-touching
vector, then classify what the surviving atoms are and how each would be set.
"""
import os, sys, collections, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
SEED = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
av = L.all_atom_values(v)
FREE = set(ad.FREE)
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
rows = sorted(OB); cols = sorted(ACTIVE); ci = {a: j for j, a in enumerate(cols)}
n, m = len(rows), len(cols)
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
settable = set()
handle = {}
for a in cols:
    for w in set(L.avars[a]):
        if w in FREE: settable.add(a); handle[a] = ('direct', w); break
    if a in settable: continue
    done = False
    for w in set(L.avars[a]):
        d = L.definer.get(w)
        if d is None: continue
        for u in set(L.avars[d]):
            if u in FREE and dz(d, u) and dz(d, u) % P == 0:
                settable.add(a); handle[a] = ('p-handle', u); done = True; break
        if done: break
non = [a for a in cols if a not in settable]
MM = []
for e in rows:
    mm, sq, co = L.eq_atoms[e]
    MM.append([co.get(a, 0) % P for a in cols])
for a in non:
    r = [0] * m; r[ci[a]] = 1; MM.append(r)
nn = len(MM); piv, r_ = [], 0
for j in range(m):
    k = next((i for i in range(r_, nn) if MM[i][j]), None)
    if k is None: continue
    MM[r_], MM[k] = MM[k], MM[r_]
    inv = pow(MM[r_][j], -1, P)
    MM[r_] = [x * inv % P for x in MM[r_]]
    for i in range(nn):
        if i != r_ and MM[i][j]:
            f = MM[i][j]
            MM[i] = [(a2 - f * b2) % P for a2, b2 in zip(MM[i], MM[r_])]
    piv.append(j); r_ += 1
ps = set(piv)
B = []
for fc in [j for j in range(m) if j not in ps]:
    z = [0] * m; z[fc] = 1
    for i, pj in enumerate(piv): z[pj] = (-MM[i][fc]) % P
    B.append(z)
print(f'kernel dim {len(B)}; supports {[sum(1 for x in b if x) for b in B]}', flush=True)
seedj = [ci[a] for a in SEED]

def supp(z): return sum(1 for x in z if x)
best = None
for start in range(len(B)):
    z = B[start][:]
    if not any(z[j] for j in seedj): continue
    improved = True
    while improved:
        improved = False
        for k in range(m):
            if not z[k]: continue
            for b in B:
                if not b[k]: continue
                lam = (-z[k]) * pow(b[k], -1, P) % P
                w = [(z[i] + lam * b[i]) % P for i in range(m)]
                if any(w[j] for j in seedj) and supp(w) < supp(z):
                    z = w; improved = True; break
            if improved: break
    if best is None or supp(z) < supp(best): best = z
sp = [cols[j] for j in range(m) if best[j]]
print(f'\nSPARSEST seed-touching kernel vector: support {len(sp)}')
print(f'  seed atoms in it: {[a for a in SEED if a in sp]}')
kinds = collections.Counter(handle.get(a, ('none',))[0] for a in sp)
print(f'  how each would be set: {dict(kinds)}')
print(f'  currently nonzero: {[a for a in sp if av[a]]}')
print(f'  equations touched: {len(set().union(*[set(L.atom2eq[a]) for a in sp]))}')
print(f'\n  atoms: {sp[:60]}')
import json
json.dump({'support': sp, 'values': [str(best[ci[a]]) for a in sp]},
          open(os.path.join(HERE, 'kervec.json'), 'w'))
print('  saved kervec.json')
