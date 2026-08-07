"""S11 step 39: the enlarged atom set SEVEN + {a22231, a22232, a22233}.

a37887 = Q^2 with Q = a22231 - 3x_18253 - 9x_23754 - ... , so moving x_28730 by d
can be undone by moving x_23754 by -d/9 -- but x_23754 is pinned by a22232, so
that atom (and its partner a22233) joins the residual.  Compute the exact optimum
for the enlarged set and compare with 7.
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
for EXTRA, cong in [([22231], 2), ([22231, 22232, 22233], 3),
                    ([22231, 22232, 22233, 22234, 22235], 4),
                    ([22231, 22232, 22233, 22234, 22235, 19087, 19088], 5)]:
    ATOMS = SEVEN + EXTRA
    Eset = set()
    for a in ATOMS: Eset |= set(L.atom2eq[a])
    E = sorted(Eset)
    rows = []
    for e in E:
        m, sq, co = L.eq_atoms[e]
        rows.append([co.get(a, 0) for a in ATOMS])
    n = len(ATOMS)
    def kdim(sel):
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
    bestk = 0
    for k in range(len(E), 0, -1):
        if k > 13: continue          # keep the enumeration finite
        hit = False
        for sel in itertools.combinations(range(len(E)), k):
            if kdim(sel) >= cong: hit = True; break
        if hit: bestk = k; break
    print(f'atoms {len(ATOMS)} (extra {EXTRA}), equations {len(E)}, '
          f'congruences {cong}')
    print(f'   max satisfiable (kernel dim >= {cong}): {bestk}  -> failing '
          f'{len(E) - bestk}   score {L.NEQ - (len(E) - bestk)}')
