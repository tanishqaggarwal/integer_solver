"""S10 step 73: the last untried lever -- what is the MINIMUM number of equations
that must break to free the wire root?

The wire-identity system M d = 0 (219 equations, 220 unknowns) has rank 217 and
d_root = 0 on its whole kernel, so freeing the root costs something.  The UNIFORM
shift costs 12 (the root pin's equations).  But is 12 minimal?

e_root lies in rowspace(M), so write  e_root = y^T M.  Any d annihilating the
support of y also has d_root = 0; hence at least one row of supp(y) must break.
The representations form  y_0 + ker(M^T), and dim ker(M^T) = 219 - 217 = 2, so we
can zero at most two coordinates.  Compute y_0, minimise its support over that
2-parameter family, and read off the true cost.

If it comes in under 7, the wire root frees for less than the give-up cost and
the trapdoor falls.
"""
import os, sys, json, collections
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T

P = 2**256 - 2**32 - 977
base = L.load(os.path.join(HERE, 'forward_state.json'))
WIRE = sorted(u for u in range(L.NVARS) if base[u] == P)
widx = {u: i for i, u in enumerate(WIRE)}
N = len(WIRE)
ROOT = widx[26064]

IDENT = {}
for a in range(L.NA):
    vs = L.avars[a]
    if not vs or not all(u in widx for u in vs):
        continue
    form = collections.defaultdict(int); ok = True
    for m, c in L.polys[a].items():
        if len(m) == 1:
            form[widx[m[0]]] += c
        elif len(m) != 0:
            ok = False; break
    if ok:
        IDENT[a] = dict(form)

EQS = sorted(set().union(*[set(L.atom2eq.get(a, ())) for a in IDENT]))
rows, rowid = [], []
for e in EQS:
    m, sq, co = L.eq_atoms[e]
    form = collections.defaultdict(int)
    for a, c in co.items():
        if a in IDENT:
            for j, cc in IDENT[a].items():
                form[j] += c * cc
    form = {j: c for j, c in form.items() if c}
    if form:
        rows.append([form.get(j, 0) for j in range(N)]); rowid.append(e)
R = len(rows)
print(f'identity system: {R} equations x {N} unknowns; root index {ROOT}')
print(f'root pin a37694 lives in {len(L.atom2eq.get(37694, {}))} equations '
      f'(the uniform-shift cost)')

# ---- solve y^T M = e_root  (i.e. M^T y = e_root) over Q --------------------
# M^T is N x R
A = [[Fraction(rows[i][j]) for i in range(R)] + [Fraction(1 if j == ROOT else 0)]
     for j in range(N)]
m, n = N, R
piv = []
r = 0
for c in range(n):
    k = next((i for i in range(r, m) if A[i][c] != 0), None)
    if k is None:
        continue
    A[r], A[k] = A[k], A[r]
    pv = A[r][c]
    A[r] = [x / pv for x in A[r]]
    for i in range(m):
        if i != r and A[i][c] != 0:
            f = A[i][c]
            A[i] = [A[i][j2] - f * A[r][j2] for j2 in range(n + 1)]
    piv.append(c); r += 1
    if r == m:
        break
inconsistent = any(all(A[i][j2] == 0 for j2 in range(n)) and A[i][n] != 0
                   for i in range(r, m))
print(f'rank(M^T) = {r}; system for y is '
      f'{"INCONSISTENT (e_root NOT in rowspace!)" if inconsistent else "consistent"}')

if inconsistent:
    print('\n*** e_root is NOT in the row space -> the kernel already moves the root '
          'for FREE. Re-check: this contradicts the measured d_root = 0.')
else:
    y = [Fraction(0)] * n
    for i, c in enumerate(piv):
        y[c] = A[i][n]
    supp = [rowid[j] for j in range(n) if y[j] != 0]
    print(f'\nparticular representation y_0: support {len(supp)} equations')
    print(f'  {supp}')
    rootqs = sorted(L.atom2eq.get(37694, ()))
    print(f'  root-pin equations: {rootqs}')
    print(f'  support == root-pin equations: {sorted(supp) == rootqs}')
    print(f'  support minus root-pin eqs: {sorted(set(supp) - set(rootqs))}')
    print(f'\n=> at least one equation in the support must break to free the root.')
    print(f'   |supp| = {len(supp)}; ker(M^T) has dim {n - r}, so at most '
          f'{n - r} coordinates can be zeroed away.')
    print(f'   MINIMUM COST >= {max(1, len(supp) - (n - r))} equations '
          f'(budget to beat: 7)')
    json.dump({'support': supp, 'rootpin_eqs': rootqs, 'kerdim': n - r},
              open(os.path.join(HERE, 'rootfree.json'), 'w'))
