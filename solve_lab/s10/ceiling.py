"""S11 step 57: the minimum-equation-cost coset leader in the WITNESS frame.

Since rank = columns in every closure measured, the solve is injective and the real
quantity is not a kernel but the coset {Mx - r}.  The achievable score is
    39033 - |equations of the atoms in supp(Mx - r)|
The canonical frame's ceiling was measured at 39,018 (below 39,026).  The witness
frame -- where the deliverable actually lives -- has never been measured.
Information-set decoding: choose |cols| rows to zero exactly, read off the cost.
"""
import os, sys, random, time, json
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd, score, grad
from fwdad import jac_column
P = ad.P
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
BAD = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(v)
vm = [x % P for x in v]
av = L.all_atom_values(v)
print(f'witness frame score {score(v)}; failing checks {BAD}', flush=True)

def jc(u):
    dv = {u: 1}
    for t in ORDER:
        a = definer[t]
        d = ad.dpart(a, t, vm)
        if d % P == 0: dv[t] = 0; continue
        s = 0
        for w in L.avars[a]:
            if w == t: continue
            dw = dv.get(w, 0)
            if dw: s += ad.dpart(a, w, vm) * dw
        dv[t] = (-s % P) * pow(d, -1, P) % P
    out = {}
    for c in CHECKS:
        s = 0
        for w in L.avars[c]:
            dw = dv.get(w, 0)
            if dw: s += ad.dpart(c, w, vm) * dw
        if s % P: out[c] = s % P
    return out

gc = {}
def gr(c):
    if c not in gc: gc[c] = set(grad(c, vm))
    return gc[c]
U, rows, cols = set(), set(BAD), {}
for it in range(9):
    nu = set()
    for c in rows: nu |= gr(c)
    nu -= U; U |= nu
    for u in sorted(nu): cols[u] = jc(u)
    nr = set(BAD)
    for u in U: nr |= set(cols[u])
    grew = nr - rows; rows = nr
    if not nu and not grew: break
rows = sorted(rows); Us = sorted(U)
print(f'witness-frame closure: {len(rows)} rows x {len(Us)} cols', flush=True)
ri = {c: i for i, c in enumerate(rows)}
n, m = len(rows), len(Us)
M = [[cols[u].get(c, 0) % P for u in Us] for c in rows]
r = [(-av[c]) % P for c in rows]

def cost_of(res):
    """Equations that ACTUALLY fail.

    Counting every equation the nonzero atoms touch overcounts badly: an equation
    holding several of them fails only if their COMBINATION is nonzero.  At the
    delivered witness the seven atoms touch 12 equations but only 7 fail.
    Evaluated mod p, so this is a lower bound on the failures over Z.
    """
    vals = {rows[i]: x for i, x in enumerate(res) if x}
    e = set()
    for a in vals:
        e |= set(L.atom2eq[a])
    bad = 0
    for eq in e:
        mm, sq, co = L.eq_atoms[eq]
        s = 0
        for a, c in co.items():
            if a in vals:
                s += c * vals[a]
        if s % P:
            bad += 1
    return bad

base = cost_of(r)
print(f'x = 0 (do nothing): residual on {sum(1 for x in r if x)} rows, '
      f'cost {base} equations -> score {L.NEQ - base}', flush=True)

best = (base, None)
t0 = time.time()
trials = int(sys.argv[1]) if len(sys.argv) > 1 else 400
budget = float(sys.argv[2]) if len(sys.argv) > 2 else 540
for t in range(trials):
    if time.time() - t0 > budget: break
    idx = random.sample(range(n), min(m, n))
    A = [M[i][:] + [r[i]] for i in idx]
    piv, r_ = [], 0
    nn = len(A)
    for j in range(m):
        k = next((i for i in range(r_, nn) if A[i][j]), None)
        if k is None: continue
        A[r_], A[k] = A[k], A[r_]
        inv = pow(A[r_][j], -1, P)
        A[r_] = [y * inv % P for y in A[r_]]
        for i in range(nn):
            if i != r_ and A[i][j]:
                f = A[i][j]
                A[i] = [(y - f * z) % P for y, z in zip(A[i], A[r_])]
        piv.append(j); r_ += 1
    if any(A[i][m] for i in range(r_, nn)): continue
    x = [0] * m
    for i, j in enumerate(piv): x[j] = A[i][m]
    res = [(sum(M[i][j] * x[j] for j in range(m)) - r[i]) % P for i in range(n)]
    c = cost_of(res)
    if c < best[0]:
        best = (c, x)
        print(f'  trial {t}: cost {c} equations -> ceiling {L.NEQ - c}', flush=True)
print(f'\nWITNESS-FRAME LINEAR CEILING: {L.NEQ - best[0]}  (cost {best[0]} equations)')
print(f'  deliverable is 39,026; canonical-frame ceiling was 39,018')
json.dump({'ceiling': L.NEQ - best[0], 'cost': best[0], 'rows': len(rows),
           'cols': len(Us)}, open(os.path.join(HERE, 'ceiling_f2.json'), 'w'))
