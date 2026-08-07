"""S11 step 89: solve ALL the remaining congruences at once, at the equation level.

Closing them one at a time oscillates -- abfix zeroes the three primitives exactly
and breaks a688/a1618; advgraph re-solves the advice and the pair comes back.  That
is the signature of a system being solved coordinate-wise when it needs to be solved
jointly.

So build the joint system at the current state:

  unknowns  every free input that structurally reaches a nonzero check, minus the
            handles (invisible mod p, and the integer lift deals with them anyway)
  rows      EVERY equation those inputs can touch -- the combination S_e, one row per
            equation, and none at all from a satisfied square (§127).  Failing
            equations get their residual as right-hand side; holding ones get 0.

If the affine system is consistent, its solution satisfies everything at once and the
lift finishes it over Z.  If it is inconsistent, the witnessing rows name exactly
which equations cannot be reconciled with the rest in one linear step -- and because
§122 says the map is a genuine polynomial, that is a statement about the tangent
space only, so the witnesses are the right targets for a finite jump.

Usage: jsolve.py [state.json]
"""
import os, sys, time
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
from intad import jacZ
import suppfree
P = ad.P
src = sys.argv[1] if len(sys.argv) > 1 else 'AB_39014.json'
v0 = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v0, rounds=6)
vm = [x % P for x in v0]
av = L.all_atom_values(v0)
BASE = L.NEQ - len(L.failing_eqs(av))
NZ = [a for a in range(L.NA) if a not in L.atom_out and av[a]]
FAIL = set(L.failing_eqs(av))
print('%s: score %d; nonzero checks %s; failing equations %d'
      % (src, BASE, NZ, len(FAIL)), flush=True)
CHECKS = [a for a in range(L.NA) if a not in L.atom_out]
_, freelist, SVS = suppfree.build(v0, modp=None)
U = set()
for c in NZ:
    m = suppfree.atom_supp(c, v0, SVS, modp=None)
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
U = sorted(U)
t0 = time.time()
cols = {u: jac_column(u, v0, vm, CHECKS) for u in U}
U = [u for u in U if cols[u]]
print('%d non-handle knobs (%.0fs)' % (len(U), time.time() - t0), flush=True)
atoms_moved = sorted(set().union(*[set(cols[u]) for u in U]))
EQS = sorted(set().union(*[set(L.atom2eq[a]) for a in atoms_moved]) | FAIL)
print('%d atoms move -> %d equations in the system' % (len(atoms_moved), len(EQS)),
      flush=True)


def comb(a_vals, e):
    s = 0
    for a, c in L.eq_atoms[e][2].items():
        if a_vals[a]:
            s += c * a_vals[a]
    return s % P


rows, rhs, rid = [], [], []
for e in EQS:
    r = []
    nzr = False
    for u in U:
        s = 0
        for a, c in L.eq_atoms[e][2].items():
            d = cols[u].get(a, 0)
            if d:
                s += c * d
        s %= P
        r.append(s)
        nzr = nzr or s
    b = (-comb(av, e)) % P
    if not nzr and not b:
        continue
    rows.append(r)
    rhs.append(b)
    rid.append(e)
print('system: %d rows x %d cols  (%d of the rows currently FAIL)'
      % (len(rows), len(U), sum(1 for e in rid if e in FAIL)), flush=True)

A = [rows[i][:] + [rhs[i]] for i in range(len(rows))]
ids = list(rid)
m, n = len(U), len(A)
piv, r_ = [], 0
for j in range(m):
    k = next((i for i in range(r_, n) if A[i][j]), None)
    if k is None:
        continue
    A[r_], A[k] = A[k], A[r_]
    ids[r_], ids[k] = ids[k], ids[r_]
    inv = pow(A[r_][j], -1, P)
    A[r_] = [x * inv % P for x in A[r_]]
    for i in range(n):
        if i != r_ and A[i][j]:
            f = A[i][j]
            A[i] = [(x - f * z) % P for x, z in zip(A[i], A[r_])]
    piv.append(j)
    r_ += 1
bad = [ids[i] for i in range(r_, n) if A[i][m]]
print('rank %d; INCONSISTENT ROWS %d' % (r_, len(bad)), flush=True)
if bad:
    print('  witnesses (equations): %s' % bad[:25])
    print('  of those, currently failing: %s' % [e for e in bad if e in FAIL][:25])
else:
    d = [0] * m
    for i, j in enumerate(piv):
        d[j] = A[i][m]
    v = list(v0)
    for j, u in enumerate(U):
        if d[j]:
            v[u] = v[u] + d[j]
    ad.fwd(v, rounds=6)
    aw = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(aw))
    print('*** CONSISTENT -- applied, score %d, nonzero checks %s'
          % (s, [a for a in CHECKS if aw[a]]), flush=True)
    _, fl2, S2 = suppfree.build(v, modp=None)
    for _ in range(8):
        aw = L.all_atom_values(v)
        todo = [a for a in CHECKS if aw[a] and aw[a] % P == 0]
        if not todo:
            break
        moved = False
        for c in todo:
            mm = suppfree.atom_supp(c, v, S2, modp=None)
            for i in range(len(fl2)):
                if not ((mm >> i) & 1):
                    continue
                u = fl2[i]
                g = jacZ(u, v, [c]).get(c, 0)
                if not g or g % P or aw[c] % g:
                    continue
                w = list(v)
                w[u] = w[u] - aw[c] // g
                ad.fwd(w, rounds=6)
                a2 = L.all_atom_values(w)
                if a2[c] == 0:
                    v, aw, moved = w, a2, True
                    break
            if moved:
                break
        if not moved:
            break
    aw = L.all_atom_values(v)
    s = L.NEQ - len(L.failing_eqs(aw))
    print('after the lift: score %d' % s)
    T.save(v, os.path.join(HERE, 'JS_%d.json' % s))
    print('saved JS_%d.json' % s)
