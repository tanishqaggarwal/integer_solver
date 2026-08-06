"""S11 step 108: coset leader with a COMPENSATION column set.

ceil15 built its closure from the free inputs that reach the FAILING CHECKS, and got
39,017.  But an equation is satisfied when its atom COMBINATION vanishes, so the way
to save an equation is not only to zero the atom that is wrong -- it is also to make
some OTHER atom of that equation nonzero in a cancelling amount.  Those compensating
atoms are moved by free inputs that never touch a failing check at all, so ceil15's
column set could not see them.

That is exactly the mechanism the 39,026 deliverable runs on (§152): seven nonzero
atoms arranged so that five of the twelve equations they touch cancel.  It was found
by search.  Here it is set up as linear algebra:

  rows     every equation touched by the enlarged column set -- the 18 that currently
           FAIL want their combination driven to zero, the rest want to stay zero
  columns  every free input that moves ANY atom of ANY currently-failing equation

then information-set decoding for the assignment leaving the fewest equations broken.

Usage: compens.py [TRIALS] [BUDGET] [state.json]
"""
import os, sys, random, time, json
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
from fwdad import jac_column
import suppfree
P = ad.P
random.seed(1013)
TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 500
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 1200
src = sys.argv[3] if len(sys.argv) > 3 else 'PF_best_39015.json'
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
vm = [x % P for x in v]
av = L.all_atom_values(v)
BASE = L.NEQ - len(L.failing_eqs(av))
FAIL = sorted(L.failing_eqs(av))
print('%s: score %d; %d failing equations' % (src, BASE, len(FAIL)), flush=True)

# columns: every free input moving ANY atom of ANY failing equation
_, freelist, SVS = suppfree.build(v, modp=None)
ATOMS = sorted({a for e in FAIL for a in L.eq_atoms[e][2]})
print('%d atoms appear in the failing equations' % len(ATOMS), flush=True)
U = set()
for a in ATOMS:
    m = suppfree.atom_supp(a, v, SVS, modp=None)
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
U = sorted(U)
print('%d free inputs move at least one of them' % len(U), flush=True)
CHECKS = [a for a in range(L.NA) if a not in L.atom_out]
t0 = time.time()
cols = {}
for i, u in enumerate(U):
    cols[u] = jac_column(u, v, vm, CHECKS)
    if i % 100 == 0:
        print('   column %d/%d (%.0fs)' % (i, len(U), time.time() - t0), flush=True)
U = [u for u in U if cols[u]]
moved = sorted(set().union(*[set(cols[u]) for u in U]))
EQS = sorted(set().union(*[set(L.atom2eq[a]) for a in moved]) | set(FAIL))
print('columns %d, atoms moved %d, equations in the system %d (%.0fs)'
      % (len(U), len(moved), len(EQS), time.time() - t0), flush=True)


def comb(a_vals, e):
    s = 0
    for a, c in L.eq_atoms[e][2].items():
        if a_vals[a]:
            s += c * a_vals[a]
    return s % P


rows, rhs = [], []
for e in EQS:
    r = []
    nz = False
    for u in U:
        s = 0
        for a, c in L.eq_atoms[e][2].items():
            d = cols[u].get(a, 0)
            if d:
                s += c * d
        s %= P
        r.append(s)
        nz = nz or s
    b = (-comb(av, e)) % P
    if not nz and not b:
        continue
    rows.append(r)
    rhs.append(b)
n, m = len(rows), len(U)
print('system: %d equation rows x %d columns' % (n, m), flush=True)


def cost(x):
    bad = 0
    for i in range(n):
        s = rhs[i]
        for j in range(m):
            if x[j] and rows[i][j]:
                s -= rows[i][j] * x[j]
        if s % P:
            bad += 1
    return bad


base = cost([0] * m)
print('do nothing: %d of these rows fail (checker says %d)' % (base, L.NEQ - BASE),
      flush=True)
best = (base, None)
t0 = time.time()
for t in range(TRIALS):
    if time.time() - t0 > BUDGET:
        break
    idx = random.sample(range(n), min(m, n))
    A = [rows[i][:] + [rhs[i]] for i in idx]
    nn = len(A)
    piv, r_ = [], 0
    for j in range(m):
        k = next((i for i in range(r_, nn) if A[i][j]), None)
        if k is None:
            continue
        A[r_], A[k] = A[k], A[r_]
        inv = pow(A[r_][j], -1, P)
        A[r_] = [y * inv % P for y in A[r_]]
        for i in range(nn):
            if i != r_ and A[i][j]:
                f = A[i][j]
                A[i] = [(y - f * z) % P for y, z in zip(A[i], A[r_])]
        piv.append(j)
        r_ += 1
    if any(A[i][m] for i in range(r_, nn)):
        continue
    x = [0] * m
    for i, j in enumerate(piv):
        x[j] = A[i][m]
    c = cost(x)
    if c < best[0]:
        best = (c, x)
        print('  trial %d: %d rows fail -> ceiling %d' % (t, c, L.NEQ - c), flush=True)
print('\nCOMPENSATION CEILING at %s : %d  (cost %d)  [%d trials, %.0fs]'
      % (src, L.NEQ - best[0], best[0], t, time.time() - t0))
print('  ceil15 (failing-check columns only) gave 39,017; the deliverable is 39,026')
json.dump({'state': src, 'ceiling': L.NEQ - best[0], 'cols': m, 'rows': n},
          open(os.path.join(HERE, 'compens.json'), 'w'))
