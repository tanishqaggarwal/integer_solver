"""S11 step 47: kernel restricted to SETTABLE atoms.

An atom is settable if it has a free variable or a p-absorbing free handle; the
rest can only be moved indirectly.  Compute  ker(M) intersect {z = 0 on the
non-settable atoms}  and ask whether it still touches the seed.  If it does, every
atom that must go nonzero is one we can actually move.
"""
import os, sys, collections, time
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
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
rows = sorted(OB); cols = sorted(ACTIVE)
ci = {a: j for j, a in enumerate(cols)}
n, m = len(rows), len(cols)
print(f'closure {n} x {m}', flush=True)

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
    ok = False
    for w in set(L.avars[a]):
        if w in FREE:
            ok = True; break
    if not ok:
        # p-absorbing handle through the atom that defines one of its variables
        for w in set(L.avars[a]):
            d = L.definer.get(w)
            if d is None: continue
            for u in set(L.avars[d]):
                if u in FREE and dz(d, u) and dz(d, u) % P == 0: ok = True; break
            if ok: break
    if ok: settable.add(a)
non = [a for a in cols if a not in settable]
print(f'settable atoms {len(settable)}; NON-settable {len(non)}')
print(f'  seed all settable? {[a for a in SEED if a not in settable]}')

# rows of M, plus one row per non-settable atom forcing it to zero
MM = []
for e in rows:
    mm, sq, co = L.eq_atoms[e]
    MM.append([co.get(a, 0) % P for a in cols])
for a in non:
    r = [0] * m; r[ci[a]] = 1
    MM.append(r)
nn = len(MM)
piv, r_ = [], 0
t0 = time.time()
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
print(f'augmented system {nn} x {m}: rank {r_}, kernel dim {m - r_}  '
      f'({time.time()-t0:.0f}s)')
ps = set(piv)
hits = []
for fc in [j for j in range(m) if j not in ps]:
    z = [0] * m; z[fc] = 1
    for i, pj in enumerate(piv): z[pj] = (-MM[i][fc]) % P
    s = [cols[j] for j in [ci[a] for a in SEED] if z[j]]
    if s: hits.append((cols[fc], s, sum(1 for x in z if x),
                       [cols[j] for j in range(m) if z[j]]))
print(f'\nkernel vectors (settable-only) touching the seed: {len(hits)}')
for fc, s, supp, full in hits[:6]:
    print(f'  free col a{fc}: seed {s}; support {supp}')
if not hits:
    print('  NONE -- forcing the non-settable atoms to zero kills every '
          'seed-touching kernel vector')
    print(f'  => the compensation MUST pass through the {len(non)} non-settable atoms')
