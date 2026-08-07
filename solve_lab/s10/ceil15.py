"""S11 step 105: the coset-leader ceiling AT the new 39,015 state.

§152 separated two different optima.  39,026 is a CODING optimum -- seven nonzero
atoms arranged so their bundle combinations cancel in all but seven equations --
while 39,015 is the ALGEBRAIC optimum, where the two-condition primitive actually closes.  The
question that decides whether the algebraic path can ever overtake the coding one is
how much coding headroom 39,015 itself has.

So run the coset leader there: build the closure of the failing checks, and ask for
the assignment of the reachable atom values that leaves the fewest EQUATIONS failing
(counting an equation as failing only when its atom COMBINATION is nonzero -- the
correction of §57).  Information-set decoding over the closure.

Usage: ceil15.py [TRIALS] [BUDGET] [state.json]
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
random.seed(97)
TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 400
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 900
src = sys.argv[3] if len(sys.argv) > 3 else 'PF_best_39015.json'
v = L.load(src if os.path.isabs(src) else os.path.join(HERE, src))
ad.fwd(v, rounds=6)
vm = [x % P for x in v]
av = L.all_atom_values(v)
BASE = L.NEQ - len(L.failing_eqs(av))
CHECKS = [a for a in range(L.NA) if a not in L.atom_out]
BAD = [a for a in CHECKS if av[a]]
print('%s: score %d; failing checks %s' % (src, BASE, BAD), flush=True)

_, freelist, SVS = suppfree.build(v, modp=None)
U = set()
for c in BAD:
    m = suppfree.atom_supp(c, v, SVS, modp=None)
    U |= {freelist[i] for i in range(len(freelist)) if (m >> i) & 1}
U = sorted(U)
t0 = time.time()
cols = {u: jac_column(u, v, vm, CHECKS) for u in U}
U = [u for u in U if cols[u]]
rows = sorted(set().union(*[set(cols[u]) for u in U]) | set(BAD))
print('closure: %d rows x %d cols (%.0fs)' % (len(rows), len(U), time.time() - t0),
      flush=True)
n, m = len(rows), len(U)
M = [[cols[u].get(c, 0) % P for u in U] for c in rows]
r = [(-av[c]) % P for c in rows]


def cost_of(res):
    """Equations whose atom COMBINATION is nonzero (§57's correction)."""
    vals = {rows[i]: x for i, x in enumerate(res) if x}
    e = set()
    for a in vals:
        e |= set(L.atom2eq[a])
    bad = 0
    for eq in e:
        s = 0
        for a, c in L.eq_atoms[eq][2].items():
            if a in vals:
                s += c * vals[a]
        if s % P:
            bad += 1
    return bad


base = cost_of(r)
print('do nothing: %d equations fail -> score %d (checker says %d)'
      % (base, L.NEQ - base, BASE), flush=True)
best = (base, None)
t0 = time.time()
for t in range(TRIALS):
    if time.time() - t0 > BUDGET:
        break
    idx = random.sample(range(n), min(m, n))
    A = [M[i][:] + [r[i]] for i in idx]
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
    res = [(sum(M[i][j] * x[j] for j in range(m)) - r[i]) % P for i in range(n)]
    c = cost_of(res)
    if c < best[0]:
        best = (c, x)
        print('  trial %d: %d equations -> ceiling %d' % (t, c, L.NEQ - c),
              flush=True)
print('\nLINEAR CEILING AT %s : %d  (cost %d equations)'
      % (src, L.NEQ - best[0], best[0]))
print('  the deliverable is 39,026; this state is %d' % BASE)
json.dump({'state': src, 'ceiling': L.NEQ - best[0], 'base': BASE,
           'rows': len(rows), 'cols': len(U)},
          open(os.path.join(HERE, 'ceil15.json'), 'w'))
