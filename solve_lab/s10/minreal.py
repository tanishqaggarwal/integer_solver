"""S11 step 53: minimise REALISABILITY cost over the settable kernel.

The sparsest vector (69 atoms) needs 42 detachments, which touch 39 outside atoms
and put 110 equations at risk.  But sparsity is the wrong objective.  Search the
8-dimensional settable kernel for the vector whose SUPPORT NEEDS THE FEWEST
DETACHMENTS / touches the fewest outside atoms.  Cost = equations of the outside
atoms not already inside the support's equations.
"""
import os, sys, collections, itertools, random, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
random.seed(97)
SEED = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
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
for a in cols:
    if any(w in FREE for w in set(L.avars[a])): settable.add(a); continue
    done = False
    for w in set(L.avars[a]):
        d0 = L.definer.get(w)
        if d0 is None: continue
        for u in set(L.avars[d0]):
            if u in FREE and dz(d0, u) and dz(d0, u) % P == 0:
                settable.add(a); done = True; break
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
seedj = [ci[a] for a in SEED]
print(f'settable kernel dimension {len(B)}', flush=True)

def cost(z):
    sp = [cols[j] for j in range(m) if z[j]]
    if not any(z[j] for j in seedj): return None
    E = set()
    for a in sp: E |= set(L.atom2eq[a])
    det = [L.atom_out[a][1] for a in sp if a in L.atom_out]
    outside = set()
    for t in det:
        for b in L.var_atoms[t]:
            if b not in set(sp): outside.add(b)
    eo = set()
    for b in outside: eo |= set(L.atom2eq[b])
    return (len(eo - E), len(sp), len(det), len(outside))

best = None
t0 = time.time()
for b in B:
    c = cost(b)
    if c and (best is None or c[0] < best[0][0]): best = (c, b)
print(f'basis vectors: best realisability cost {best[0]}')
tried = 0
for it in range(4000):
    k = random.randint(2, min(4, len(B)))
    idx = random.sample(range(len(B)), k)
    z = [0] * m
    for i in idx:
        lam = random.randrange(1, 40)
        z = [(z[j] + lam * B[i][j]) % P for j in range(m)]
    c = cost(z)
    tried += 1
    if c and c[0] < best[0][0]:
        best = (c, z); print(f'  it{it}: cost {c}', flush=True)
c, z = best
sp = [cols[j] for j in range(m) if z[j]]
print(f'\ntried {tried} combinations ({time.time()-t0:.0f}s)')
print(f'BEST: net equations at risk {c[0]}; support {c[1]} atoms; '
      f'{c[2]} detachments touching {c[3]} outside atoms')
print(f'  seed atoms in it: {[a for a in SEED if a in sp]}')
