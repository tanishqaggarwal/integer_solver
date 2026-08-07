"""S11 step 58: equation-space compensation for the 39,017 state's blockers.

jm_azero00_39017.json has only THREE nonzero atoms -- a688, a1618, a40608 -- costing
16 equations, with all twelve gadget equations satisfied.  The joint-move sweep
tried single VARIABLE moves from it.  Apply the compensation machinery instead:
which atoms have their entire equation set inside those 16, and does adding them
create a kernel?
"""
import os, sys, itertools
from fractions import Fraction
HERE='/home/user/integer_solver/solve_lab/s10'; LAB='/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB,'s9','eff')); sys.path.insert(0, HERE)
import lib as L
S = [688, 1618, 40608]
E = sorted(set().union(*[set(L.atom2eq[a]) for a in S]))
ES = set(E)
print(f'blockers {S} touch {len(E)} equations: {E}')
for a in S: print(f'  a{a}: {len(L.atom2eq[a])} equations')
inside = {}
for e in E:
    m, sq, co = L.eq_atoms[e]
    for a, c in co.items():
        if c: inside[a] = inside.get(a, 0) + 1
cands = []
for a, k in inside.items():
    if a in S: continue
    out = len(set(L.atom2eq[a]) - ES)
    cands.append((out, k, a))
cands.sort()
print(f'\natoms appearing in those equations: {len(inside)}')
print(f'compensator candidates (outside-equations, in-how-many, atom):')
for out, k, a in cands[:14]:
    print(f'  a{a:<6} in {k:>2} of the {len(E)}; {out:>3} equations OUTSIDE')
free_comp = [a for out, k, a in cands if out == 0]
print(f'\ncompensators with ZERO equations outside: {free_comp}')

def rank_q(cols_):
    if not cols_: return 0
    n = len(cols_[0]); m = len(cols_)
    M = [[Fraction(cols_[j][i]) for j in range(m)] for i in range(n)]
    r_ = 0
    for j in range(m):
        k = next((i for i in range(r_, n) if M[i][j] != 0), None)
        if k is None: continue
        M[r_], M[k] = M[k], M[r_]
        pv = M[r_][j]; M[r_] = [x / pv for x in M[r_]]
        for i in range(n):
            if i != r_ and M[i][j] != 0:
                f = M[i][j]
                M[i] = [x - f * y for x, y in zip(M[i], M[r_])]
        r_ += 1
    return r_
def col(a): return [L.eq_atoms[e][2].get(a, 0) for e in E]
base = [col(a) for a in S]
print(f'\nrank of the {len(E)} x {len(S)} blocker matrix: {rank_q(base)}')
pool = [a for out, k, a in cands if out <= 4][:16]
print(f'searching subsets of the {len(pool)} cheapest compensators')
best = None
for k in range(1, 6):
    for T in itertools.combinations(pool, k):
        eqs = set()
        for a in T: eqs |= (set(L.atom2eq[a]) - ES)
        if best and len(eqs) >= best[0]: continue
        cols_ = base + [col(a) for a in T]
        if rank_q(cols_) < len(S) + k:
            if best is None or len(eqs) < best[0]:
                best = (len(eqs), T, k)
                print(f'  k={k}: {T} -> kernel; {len(eqs)} equations outside')
print(f'\nBEST: {best}')
