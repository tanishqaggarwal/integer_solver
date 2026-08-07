"""S11 step 59: full compensation analysis for the 39,017 blockers.

16 equations, 25 atoms appearing in them -> a kernel of dimension >= 9 must exist.
The question is only the outside cost of the atoms it needs.  Minimise
|union of outside equations| over kernel vectors touching a688/a1618/a40608.
"""
import os, sys, itertools, random, math
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
random.seed(5)
S = [688, 1618, 40608]
E = sorted(set().union(*[set(L.atom2eq[a]) for a in S])); ES = set(E)
atoms = sorted({a for e in E for a, c in L.eq_atoms[e][2].items() if c})
ci = {a: j for j, a in enumerate(atoms)}
n, m = len(E), len(atoms)
print(f'system {n} equations x {m} atoms; blockers at columns '
      f'{[ci[a] for a in S]}')
M = [[Fraction(L.eq_atoms[e][2].get(a, 0)) for a in atoms] for e in E]
piv, r_ = [], 0
for j in range(m):
    k = next((i for i in range(r_, n) if M[i][j] != 0), None)
    if k is None: continue
    M[r_], M[k] = M[k], M[r_]
    pv = M[r_][j]; M[r_] = [x / pv for x in M[r_]]
    for i in range(n):
        if i != r_ and M[i][j] != 0:
            f = M[i][j]
            M[i] = [x - f * y for x, y in zip(M[i], M[r_])]
    piv.append(j); r_ += 1
print(f'rank {r_}, kernel dimension {m - r_}')
ps = set(piv)
B = []
for fc in [j for j in range(m) if j not in ps]:
    z = [Fraction(0)] * m; z[fc] = Fraction(1)
    for i, pj in enumerate(piv): z[pj] = -M[i][fc]
    B.append(z)
out_eqs = {a: set(L.atom2eq[a]) - ES for a in atoms}
sj = [ci[a] for a in S]
def cost(z):
    if not any(z[j] != 0 for j in sj): return None
    e = set()
    for j in range(m):
        if z[j] != 0: e |= out_eqs[atoms[j]]
    return len(e)
best = None
for b in B:
    c = cost(b)
    if c is not None and (best is None or c < best[0]): best = (c, b)
print(f'basis vectors touching the blockers: '
      f'{sum(1 for b in B if cost(b) is not None)} of {len(B)}; best cost {best[0] if best else None}')
for it in range(20000):
    k = random.randint(2, min(5, len(B)))
    z = [Fraction(0)] * m
    for i in random.sample(range(len(B)), k):
        lam = Fraction(random.randint(-6, 6))
        if lam == 0: continue
        z = [z[j] + lam * B[i][j] for j in range(m)]
    c = cost(z)
    if c is not None and (best is None or c < best[0]):
        best = (c, z); print(f'  it{it}: cost {c}')
c, z = best
sp = [atoms[j] for j in range(m) if z[j] != 0]
print(f'\nBEST kernel vector: {c} equations outside, support {len(sp)} atoms')
print(f'  support: {sp}')
print(f'  blockers in it: {[a for a in S if a in sp]}')
print(f'  => a full repair of this state would cost {c} equations '
      f'-> score {L.NEQ - c}   (deliverable 39026)')
