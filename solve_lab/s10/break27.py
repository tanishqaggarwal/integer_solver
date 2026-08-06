"""S11 step 18: A1 = 0 (a22230 = 0), repair a7930, realise A2..A6 optimally.

eq 29125 = A1 alone, so A1 = 0 buys a sixth satisfied equation -> 6 failing.
Cost of A1 = 0 is that x_28730 moves, breaking a7930 (+ its 1-eq shadow a41512).
a7930 = 9367949*(x_24548 - x_25442) - x_7927 and x_24548 is a FREE input, so its
congruence is repairable.  Then choose A2..A6 (all free, zero collateral) to kill
five of the remaining ten equations.
"""
import os, sys, itertools
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
from frame2 import definer, ORDER, FREE, CHECKS, fwd, score, grad
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
SSET = set(SEVEN)
base = L.load(os.path.join(LAB, 'best', 'new_instance_partial_39026.json'))
fwd(base)
print(f'delivered witness {score(base)}')

v = list(base)
v[28730] = v[9413] * P                     # A1 = 0
fwd(v, rounds=8)
av = L.all_atom_values(v)
print(f'after A1=0: nonzero {[a for a in range(L.NA) if av[a]]}  score {score(v)}')

# --- repair a7930 through the free input x_24548 -----------------------------
print(f'\nx_24548 consumers: {sorted(L.var_atoms[24548])}')
tgt = T.solve_lin(7930, 24548, v)
print(f'  solve a7930 for x_24548 -> {"OK" if tgt is not None else "not divisible"}')
if tgt is None:
    vm = [x % P for x in v]
    g = grad(7930, vm)
    d = g.get(24548, 0)
    r = av[7930] % P
    if d % P:
        v[24548] = v[24548] + (-r * pow(d, -1, P)) % P
        fwd(v, rounds=8)
        av = L.all_atom_values(v)
        print(f'  Newton on x_24548: a7930 mod p now {av[7930] % P == 0}')
        tgt = T.solve_lin(7930, 24548, v)
        print(f'  re-solve -> {"OK" if tgt is not None else "still not divisible"}')
if tgt is not None:
    v[24548] = tgt
    fwd(v, rounds=8)
av = L.all_atom_values(v)
nz = [a for a in range(L.NA) if av[a]]
print(f'after a7930 repair: nonzero {nz}  score {score(v)}')
outside = [a for a in nz if a not in SSET]
print(f'  atoms outside the seven: {outside}')

# --- now choose A2..A6 to kill five of the remaining equations ---------------
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SEVEN]))
rows = []
for e in E:
    m, sq, co = L.eq_atoms[e]
    rows.append([co.get(a, 0) for a in SEVEN])
A0 = av[22229]
print(f'\nA0 = {str(A0)[:26]}  A1 = {av[22230]}')

def solve_free(sel):
    """choose A2..A6 (A0 fixed, A1 = 0) killing the equations in sel; None if impossible."""
    M, b = [], []
    for i in sel:
        r = rows[i]
        M.append([Fraction(r[k]) for k in range(2, 7)])
        b.append(Fraction(-r[0] * A0))
    n, m = len(M), 5
    aug = [M[i] + [b[i]] for i in range(n)]
    piv, r_ = [], 0
    for j in range(m):
        k = next((i for i in range(r_, n) if aug[i][j] != 0), None)
        if k is None: continue
        aug[r_], aug[k] = aug[k], aug[r_]
        pv = aug[r_][j]; aug[r_] = [x / pv for x in aug[r_]]
        for i in range(n):
            if i != r_ and aug[i][j] != 0:
                f = aug[i][j]
                aug[i] = [x - f * y for x, y in zip(aug[i], aug[r_])]
        piv.append(j); r_ += 1
    for i in range(r_, n):
        if aug[i][m] != 0: return None
    sol = [Fraction(0)] * 5
    for i, j in enumerate(piv): sol[j] = aug[i][m]
    if any(x.denominator != 1 for x in sol): return None
    return [int(x) for x in sol]

best = None
for k in (6, 5):
    for sel in itertools.combinations(range(len(E)), k):
        s = solve_free(sel)
        if s is None: continue
        best = (k, sel, s); break
    if best: break
print(f'best killable set: {best[0] if best else 0} equations '
      f'{[E[i] for i in best[1]] if best else []}')
if best:
    print(f'  A2..A6 = {[str(x)[:20] for x in best[2]]}')
