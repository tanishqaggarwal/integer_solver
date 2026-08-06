"""S11 step 41: kernel DIMENSION, not just nonzero-ness.

The balance law assumed k = n - c, i.e. that the k rows are independent.  A
DEPENDENT subset has a bigger kernel, so a 7-subset of rank 6 would give kernel
dim 2 -- enough for BOTH congruences -- and 7 of 12 satisfied.
Check every subset size against the required dimension.
"""
import os, sys, itertools, collections
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
E = sorted(set().union(*[set(L.atom2eq[a]) for a in SEVEN]))

def kdim(rows, sel, n):
    M = [[Fraction(rows[i][k]) for k in range(n)] for i in sel]
    nn = len(M); r_ = 0
    for j in range(n):
        k = next((i for i in range(r_, nn) if M[i][j] != 0), None)
        if k is None: continue
        M[r_], M[k] = M[k], M[r_]
        pv = M[r_][j]; M[r_] = [x / pv for x in M[r_]]
        for i in range(nn):
            if i != r_ and M[i][j] != 0:
                f = M[i][j]
                M[i] = [x - f * y for x, y in zip(M[i], M[r_])]
        r_ += 1
    return n - r_

for ATOMS, cong, tag in [(SEVEN, 2, 'the seven'),
                         (SEVEN + [22231], 2, '+ a22231')]:
    n = len(ATOMS)
    rows = []
    for e in E:
        m, sq, co = L.eq_atoms[e]
        rows.append([co.get(a, 0) for a in ATOMS])
    print(f'\n=== {tag}: n = {n}, congruences = {cong} ===')
    for k in range(12, 2, -1):
        dist = collections.Counter()
        bestsel = {}
        for sel in itertools.combinations(range(12), k):
            d = kdim(rows, sel, n)
            dist[d] += 1
            if d not in bestsel: bestsel[d] = sel
        mx = max(dist)
        print(f'  size {k:>2}: kernel-dim distribution {dict(sorted(dist.items()))}'
              f'   max {mx}')
        if mx >= cong:
            print(f'    *** size {k} reaches kernel dim {mx} >= {cong}: '
                  f'eqs {[E[i] for i in bestsel[mx]]}')
            print(f'    -> {12-k} of the twelve fail'
                  f'{" (+1 for a37887)" if n > 7 else ""}')
            break
