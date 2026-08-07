"""How many disjoint rank-20 subsets can the 162 redundant rows be packed into?

t disjoint such subsets => any B of < t redundant rows leaves one intact => rank(A_SAT\B)
unchanged => the integer feasible set unchanged => breaking those rows is worthless.
Also measures WHY: the rank profile of R162 and the coordinate structure that limits it.
"""
import sys, os, random
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S

KN = S.KNOB
ESS = [2554, 6816, 8124, 9123, 9421, 'S']
R162 = [e for e in S.SAT if e not in ESS]
V = {e: [Fraction(S.rows[e].get(u, 0)) for u in KN] for e in S.SAT}


def basis_add(basis, w):
    w = w[:]
    for bp, bv in basis:
        if w[bp] != 0:
            f = w[bp]; w = [a - f * b for a, b in zip(w, bv)]
    p = next((i for i, x in enumerate(w) if x != 0), None)
    if p is None: return False
    d = w[p]; basis.append((p, [x / d for x in w])); return True


def rank(rows):
    b = []
    for e in rows: basis_add(b, V[e])
    return len(b)


print('rank(R162) =', rank(R162))
# greedy multi-base packing: round robin
best = 0; bestpacks = None
for trial in range(40):
    order = list(R162); random.Random(trial).shuffle(order)
    bases = []; assigned = []
    for e in order:
        placed = False
        for i, b in enumerate(bases):
            if len(assigned[i]) < 20 and basis_add(b, V[e]):
                assigned[i].append(e); placed = True; break
        if not placed:
            bases.append([]); assigned.append([])
            if basis_add(bases[-1], V[e]): assigned[-1].append(e)
    full = [a for a, b in zip(assigned, bases) if len(b) == 20]
    if len(full) > best:
        best = len(full); bestpacks = full
print('best round-robin packing: t = %d disjoint rank-20 subsets (of 162 rows, 8 is the ceiling)' % best)

# WHY: how many rows lie in a common hyperplane?  measure via coordinate sparsity in a basis of the span
b = []
basisrows = []
for e in R162:
    if basis_add(b, V[e]): basisrows.append(e)
print('a basis of span(R162) is %d rows' % len(basisrows))
# express every row in that basis -> 162 x 20 coordinate matrix; count zeros per coordinate
import itertools
B = [(p, bv) for p, bv in b]
coords = {}
for e in R162:
    w = V[e][:]; c = []
    for bp, bv in B:
        f = w[bp]; c.append(f)
        if f: w = [a - f * x for a, x in zip(w, bv)]
    coords[e] = c
for j in range(20):
    nz = sum(1 for e in R162 if coords[e][j] != 0)
    print('  coord %2d : %3d rows nonzero  -> a cocircuit of size %d if this coord is a valid functional' % (j, nz, nz))
