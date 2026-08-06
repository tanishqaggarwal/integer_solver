"""S11 step 68: AUDIT of step 66's "the exact subsystem is inconsistent".

Step 66 built a 220 x 80 exactly-linear system and found 7 inconsistent rows.
Before that is allowed to mean anything, three objections have to be answered --
the same skepticism that killed the cluster agent's infeasibility claim:

  (1) COLUMN CLOSURE.  The 80 columns were the free inputs reaching the two cluster
      targets.  The 220 rows are every check those 80 columns touch -- but a row can
      depend on free inputs OUTSIDE the 80.  Then `Mx = b` restricted to 80 columns
      is artificially over-constrained and its inconsistency says nothing.
  (2) SAMPLING.  "Exactly linear" was decided from 2 random probes per input.  Two
      agreements are weak evidence for an identity.
  (3) ATOMS vs EQUATIONS.  `Mx = b` demands every row atom be ZERO.  But an equation
      is satisfied iff its atom COMBINATION vanishes; at the delivered witness many
      atoms are nonzero while their equations hold.  Even a true inconsistency here
      would not be an infeasibility proof unless the rows sit alone in equations.

This script tests all three.  Usage: exactlin2.py [NPROBE]
"""
import os, sys, random, time, collections
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from fwdad import jac_column
import suppfree
P = ad.P
random.seed(11)
NT = int(sys.argv[1]) if len(sys.argv) > 1 else 4

v0 = L.load(os.path.join(HERE, 'mod9118_0.json'))
vm0 = [x % P for x in v0]
av0 = L.all_atom_values(v0)
CHECKS = sorted(a for a in range(L.NA) if a not in L.atom_out)
BAD = [21617, 29539]
FORBID = {2081, 4287}
U = sorted((set(ad.grad(BAD[0], vm0)) | set(ad.grad(BAD[1], vm0))) - FORBID)
print(f'base state score {L.NEQ - len(L.failing_eqs(av0))}; {len(U)} columns', flush=True)

t0 = time.time()
cols = {u: jac_column(u, v0, vm0, CHECKS) for u in U}
print(f'jacobian columns: {time.time()-t0:.0f}s', flush=True)

lin_bad = collections.Counter()
touched = set()
for u in U:
    touched |= set(cols[u])
for i, u in enumerate(U):
    col = cols[u]
    for _ in range(NT):
        d = random.randrange(1, P)
        w = list(v0)
        w[u] = w[u] + d
        ad.fwd(w, rounds=6)
        aw = L.all_atom_values(w)
        for c in col:
            if (av0[c] + col[c] * d) % P != aw[c] % P:
                lin_bad[c] += 1
    if i % 20 == 0:
        print(f'  probe {i}/{len(U)}  ({time.time()-t0:.0f}s)', flush=True)
exact = sorted(c for c in touched if lin_bad[c] == 0)
print(f'\n[2] SAMPLING: {NT} probes/input -> {len(exact)} of {len(touched)} rows '
      f'exactly linear (step 66 with 2 probes said 220)', flush=True)

idx, freelist, vs = suppfree.build(vm0)
US = set(U)
missing = {}
for c in exact:
    m = suppfree.atom_supp(c, vm0, vs)
    extra = [freelist[i] for i in range(len(freelist)) if (m >> i) & 1
             and freelist[i] not in US]
    if extra:
        missing[c] = extra
allextra = sorted(set().union(*missing.values())) if missing else []
print(f'[1] COLUMN CLOSURE: {len(missing)} of {len(exact)} rows depend on free '
      f'inputs outside the 80; {len(allextra)} distinct missing columns', flush=True)

# which rows were the inconsistent ones, and are THEY column-closed?
def solve(rows, colset):
    cl = sorted(colset)
    n, m = len(rows), len(cl)
    A = [[cols[u].get(c, 0) % P for u in cl] + [(-av0[c]) % P] for c in rows]
    piv, r_ = [], 0
    for j in range(m):
        k = next((i for i in range(r_, n) if A[i][j]), None)
        if k is None:
            continue
        A[r_], A[k] = A[k], A[r_]
        inv = pow(A[r_][j], -1, P)
        A[r_] = [x * inv % P for x in A[r_]]
        for i in range(n):
            if i != r_ and A[i][j]:
                f = A[i][j]
                A[i] = [(x - f * z) % P for x, z in zip(A[i], A[r_])]
        piv.append(j)
        r_ += 1
    return r_, [i for i in range(r_, n) if A[i][m]], A, piv, cl

r_, bad, A, piv, cl = solve(exact, U)
print(f'\nEXACT subsystem {len(exact)} x {len(U)}: rank {r_}, inconsistent {len(bad)}',
      flush=True)

closed = [c for c in exact if c not in missing]
print(f'\n--- restricted to the COLUMN-CLOSED rows only ---')
r2, bad2, _, _, _ = solve(closed, U)
print(f'{len(closed)} rows x {len(U)} cols: rank {r2}, inconsistent {len(bad2)}',
      flush=True)
if not bad2:
    print('*** step 66 was an ARTIFACT: every inconsistency came from rows whose')
    print('    columns were missing.  No infeasibility conclusion survives.')
else:
    print('*** the inconsistency SURVIVES column closure on these rows.')
    print('[3] ATOMS vs EQUATIONS: do those rows sit alone in an equation?')

# [3] equations of the exact rows: how many are size-1 (atom forced to zero)?
alone = 0
for c in closed:
    for e in L.atom2eq[c]:
        if len(L.eq_atoms[e][2]) == 1:
            alone += 1
            break
print(f'    of {len(closed)} column-closed exact rows, {alone} sit alone in some '
      f'equation (so must truly be zero)')
