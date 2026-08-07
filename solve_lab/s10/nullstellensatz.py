"""S11 step 63: search for a LINEAR INFEASIBILITY CERTIFICATE.

Every equation is literally poly_e(x) = 0, so for any solution x and any rational
y, sum_e y_e * poly_e(x) = 0.  If some combination is identically a NONZERO
CONSTANT, no solution exists -- a complete infeasibility proof, with no atom-level
or linearisation reasoning anywhere.

Method: rows = distinct monomials, columns = equations.  Find y killing every
NON-CONSTANT monomial; then look at the constant row.  Constants really do enter
(e.g. a37694 = x_26064 - p carries -p), so this is not vacuous.

Usage: nullstellensatz.py <region>   where region selects the equation set.
"""
import os, sys, collections, time
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]

def eq_poly(e):
    """the equation's polynomial: m*(S) or m*(S)^2 with S = sum c_a * atom_a"""
    m, sq, co = L.eq_atoms[e]
    S = collections.defaultdict(int)
    for a, c in co.items():
        for mono, cc in L.polys[a].items():
            S[tuple(sorted(mono))] += c * cc
    S = {k: c for k, c in S.items() if c}
    if not sq:
        return {k: m * c for k, c in S.items()}
    out = collections.defaultdict(int)
    it = list(S.items())
    for k1, c1 in it:
        for k2, c2 in it:
            out[tuple(sorted(k1 + k2))] += m * c1 * c2
    return {k: c for k, c in out.items() if c}

region = sys.argv[1] if len(sys.argv) > 1 else 'twelve'
if region == 'twelve':
    E = sorted(set().union(*[set(L.atom2eq[a]) for a in SEVEN]))
elif region == 'failing':
    v = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
    E = sorted(L.failing_eqs(L.all_atom_values(v)))
elif region == 'wire':
    E = sorted(L.atom2eq[37694])
else:
    E = sorted(set().union(*[set(L.atom2eq[int(x)]) for x in region.split(',')]))
print(f'region "{region}": {len(E)} equations {E[:14]}', flush=True)

t0 = time.time()
polys = {e: eq_poly(e) for e in E}
monos = sorted(set().union(*[set(p) for p in polys.values()]), key=lambda k: (len(k), k))
print(f'distinct monomials: {len(monos)}  ({time.time()-t0:.0f}s)', flush=True)
mi = {k: i for i, k in enumerate(monos)}
const_row = mi.get((), None)
print(f'constant monomial present: {const_row is not None}')
n, m = len(monos), len(E)
A = [[Fraction(polys[e].get(k, 0)) for e in E] for k in monos]
# eliminate over the NON-constant rows only
rows = [i for i in range(n) if monos[i] != ()]
piv, r_ = [], 0
Aw = [A[i][:] + [Fraction(1 if j == i else 0) for j in range(0)] for i in rows]
# track combination: augment with identity over columns instead (null space of A_rows)
B = [row[:] for row in [A[i] for i in rows]]
nn = len(B)
piv, r_ = [], 0
for j in range(m):
    k = next((i for i in range(r_, nn) if B[i][j] != 0), None)
    if k is None: continue
    B[r_], B[k] = B[k], B[r_]
    pv = B[r_][j]; B[r_] = [x / pv for x in B[r_]]
    for i in range(nn):
        if i != r_ and B[i][j] != 0:
            f = B[i][j]
            B[i] = [x - f * z for x, z in zip(B[i], B[r_])]
    piv.append(j); r_ += 1
print(f'rank over non-constant monomials: {r_} of {m} equations '
      f'-> null dim {m - r_}  ({time.time()-t0:.0f}s)')
ps = set(piv)
found = 0
for fc in [j for j in range(m) if j not in ps]:
    y = [Fraction(0)] * m; y[fc] = Fraction(1)
    for i, pj in enumerate(piv): y[pj] = -B[i][fc]
    c = sum(y[j] * Fraction(polys[E[j]].get((), 0)) for j in range(m))
    found += 1
    if c != 0:
        print(f'\n*** INFEASIBILITY CERTIFICATE: combination of {sum(1 for x in y if x)} '
              f'equations is identically the nonzero constant {c}')
        print(f'    equations: {[E[j] for j in range(m) if y[j]]}')
        break
else:
    if found:
        print(f'{found} null combinations found; every one has constant term 0 '
              f'(so each is the identically-zero polynomial -- no certificate here)')
    else:
        print('no null combination in this region (rank is full) -- no certificate')
