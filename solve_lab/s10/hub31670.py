"""S10 step 61: exact equation-space optimum around hub a31670.

a31670 hits FOUR of the six inconsistency certificates, so if its true cost is
small the whole hitting set is cheap.  Structure:

    a31670 = (x_22152 - HUGE) - 7550763*x_29309      [x_24601 = 1]
    a31669 = x_29309 - p*x_105                        [x_105 a free solo handle]

x_29309 and x_105 are free, so (a31670, a31669) = (D - 7550763*s, s - p*h):
a31669 is ANY integer and a31670 is fixed modulo 7550763.  But a31669 is ONE
value shared by all ten equations -- the earlier "9 of 10 compensable" count was
existence of a helper, not independence.  Solve the real system exactly.
"""
import os, sys, collections, json, math, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L
import tools as T
import ad

P = ad.P
atom_out = L.atom_out
v = L.load(os.path.join(HERE, 'forward_state.json'))
av = L.all_atom_values(v)

X = 31670
print(f'a{X}: {L.atom_src[X]}')
print(f'a31669: {L.atom_src[31669]}')

SOLO = collections.defaultdict(list)
for u in range(L.NVARS):
    if u not in L.definer and len(L.var_atoms[u]) == 1:
        SOLO[L.var_atoms[u][0]].append(u)
GRAN = {}
for a, us in SOLO.items():
    g = 0
    for u in us:
        r = T.lin_parts(a, u, v)
        if r:
            g = math.gcd(g, abs(r[0]))
    if g:
        GRAN[a] = g

EX = sorted(L.atom2eq.get(X, ()))
# candidate adjustable atoms appearing in X's equations
cand = set()
for e in EX:
    for a in L.eq_atoms[e][2]:
        if a != X and a in GRAN:
            cand.add(a)
cand = sorted(cand)
print(f'\nX equations: {len(EX)}; adjustable atoms in them: {len(cand)}')

# full equation set touched if we use those atoms
EQS = set(EX)
for a in cand:
    EQS |= set(L.atom2eq.get(a, ()))
EQS = sorted(EQS)
ATOMS = [X] + cand
print(f'touched equations if all are used: {len(EQS)}')
print(f'granularities: X free mod 7550763; helpers: '
      f'{[(a, GRAN[a]//P if GRAN[a] % P == 0 else GRAN[a]) for a in cand[:8]]} '
      f'(shown /p where divisible)')

# Matrix of the region
M = [[L.eq_atoms[e][2].get(a, 0) for a in ATOMS] for e in EQS]
n = len(ATOMS)
print(f'system {len(M)} equations x {n} atoms')

from fractions import Fraction
def rank_q(mat):
    m = [[Fraction(x) for x in r] for r in mat]
    rows, cols = len(m), len(m[0]); r = 0
    for c in range(cols):
        k = next((i for i in range(r, rows) if m[i][c] != 0), None)
        if k is None: continue
        m[r], m[k] = m[k], m[r]
        pv = m[r][c]; m[r] = [x/pv for x in m[r]]
        for i in range(rows):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [m[i][j]-f*m[r][j] for j in range(cols)]
        r += 1
        if r == rows: break
    return r

print(f'rank of the region matrix = {rank_q(M)} of {n} atom columns')
print(f'  -> with all {n} atom values free, at most {len(EQS)} equations could vanish;')
print(f'     the kernel of the transposed system bounds how many must fail.')

# how many equations can be satisfied: those in the row space sense.
# A vector of atom values y satisfies equation e iff M_e . y = 0.
# Max satisfied = |EQS| - (minimum number of independent rows violated).
# With y unconstrained in Q^n, we can satisfy any set S with rank(M_S) < n...
# but y = 0 satisfies all, so the real question is: y must have y_X != 0
# (X is being sacrificed, i.e. its value is FORCED nonzero by the chain).
print('\nIf a31670 must be nonzero, how many of the touched equations can still vanish?')
best = 0
for size in range(len(EQS), 0, -1):
    ok = False
    for S in itertools.combinations(range(len(EQS)), size):
        sub = [M[i] for i in S]
        # need a nonzero-in-coord-0 vector in the kernel of sub
        r_all = rank_q(sub)
        r_drop = rank_q([[row[j] for j in range(1, n)] for row in sub])
        if r_drop < r_all:      # X's column is independent -> y_X can be nonzero
            ok = True
            break
        if r_all < n and r_drop == r_all:
            # kernel exists but may force y_X = 0; test by rank comparison
            ok = True
            break
    if ok:
        best = size
        print(f'  max simultaneously satisfiable (with a{X} != 0): {size} of {len(EQS)}')
        print(f'  -> minimum failing in this region: {len(EQS)-size}')
        break
    if size < len(EQS) - 25:
        break
