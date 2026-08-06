"""S11 step 23: EQUATION-SPACE COMPENSATION.

The 12 equations have rank 7 over the seven residual atoms, so all 12 hold only if
A = 0.  But an extra atom b that also appears in those equations turns the matrix
into 12 x 8, which CAN have a kernel -- letting all 12 hold with A != 0.  The cost
is whatever b's OTHER equations do.  So: which atoms live inside those 12
equations, and how much of them lies outside?
"""
import os, sys, collections
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SEVEN]))
print(f'the twelve equations: {E}')
ESET = set(E)
inside = collections.Counter()
for e in E:
    m, sq, co = L.eq_atoms[e]
    for a in co: inside[a] += 1
print(f'distinct atoms appearing in them: {len(inside)}')
cands = []
for a, k in inside.items():
    if a in SEVEN: continue
    out = len(set(L.atom2eq[a]) - ESET)
    cands.append((out, k, a))
cands.sort()
print(f'\ncompensator candidates (outside-equations, in-how-many-of-the-12, atom):')
for out, k, a in cands[:25]:
    print(f'  a{a:<6} appears in {k} of the 12; {out} equations OUTSIDE')
best_out = cands[0][0] if cands else None
print(f'\ncheapest compensator lives in {best_out} equations outside the twelve')

def rank_q(rows):
    if not rows: return 0
    M = [[Fraction(x) for x in r] for r in rows]
    n, m = len(M), len(M[0]); r_ = 0
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

base_rows = []
for e in E:
    m, sq, co = L.eq_atoms[e]
    base_rows.append([co.get(a, 0) for a in SEVEN])
print(f'rank of the 12 x 7 matrix: {rank_q(base_rows)} (so A = 0 is forced)')

print('\nadding one compensator at a time -> does a kernel appear?')
hits = []
for out, k, a in cands[:40]:
    rows = []
    for e, br in zip(E, base_rows):
        m, sq, co = L.eq_atoms[e]
        rows.append(br + [co.get(a, 0)])
    r = rank_q(rows)
    if r < 8:
        hits.append((out, a, r))
        print(f'  a{a:<6} (outside {out:>3}) -> rank {r} of 8  *** KERNEL, dim {8-r}')
if not hits:
    print('  none of the cheapest 40 creates a kernel')
