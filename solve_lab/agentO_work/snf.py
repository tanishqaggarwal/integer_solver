"""Region R0+{a23618}: over Q all 13 equations are satisfiable, over Z only 6.
   Locate the integrality obstruction exactly (Smith normal form / invariant factors)."""
import sys, json, time, itertools
from fractions import Fraction
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentO_work')
import regiongrow as G, harness as H

OD = '/home/user/integer_solver/solve_lab/agentO_work'
LOG = open(OD + '/runs/snf.log', 'w', buffering=1)


def say(*a):
    print(*a, file=LOG)


R = G.R0 + [23618]
P = G.private_vars(R)
const, cols = G.build_model(R, P, G.V0)
Eqs, rows = G.eq_system(R, P, const, cols)
say('P =', P)
say('Eqs =', Eqs)

A = [[rows[e][0].get(u, 0) for u in P] for e in Eqs]
b = [rows[e][1] for e in Eqs]


def hnf_solve(A, b):
    """Row-reduce [A|b] over Z by integer row ops; report the first inconsistent pivot."""
    import copy
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    n = len(M)
    m = len(A[0])
    r = 0
    piv = []
    for c in range(m):
        rows_nz = [i for i in range(r, n) if M[i][c]]
        if not rows_nz:
            continue
        while True:
            rows_nz = [i for i in range(r, n) if M[i][c]]
            if len(rows_nz) <= 1:
                break
            rows_nz.sort(key=lambda i: abs(M[i][c]))
            p = rows_nz[0]
            for i in rows_nz[1:]:
                q = M[i][c] // M[p][c]
                for k in range(c, m + 1):
                    M[i][k] -= q * M[p][k]
        rows_nz = [i for i in range(r, n) if M[i][c]]
        if not rows_nz:
            continue
        p = rows_nz[0]
        M[r], M[p] = M[p], M[r]
        piv.append((c, M[r][c]))
        r += 1
    return M, r, piv


M, r, piv = hnf_solve(A, b)
say('rank over Z (echelon) =', r)
say('pivots (column, pivot value):')
for c, pv in piv:
    say('   x_%d : pivot %s (%d bits)' % (P[c], pv if abs(pv) < 10 ** 12 else str(pv)[:24] + '..', abs(pv).bit_length()))
say('zero rows / consistency:')
for i in range(len(M)):
    if all(x == 0 for x in M[i][:-1]):
        say('   row %d: 0 = %s  -> %s' % (i, str(M[i][-1])[:40] + ('..' if abs(M[i][-1]) > 10 ** 20 else ''),
                                          'OK' if M[i][-1] == 0 else 'INCONSISTENT over Z and Q'))
# back-substitution divisibility test
say('\nback substitution:')
sol = {}
bad = []
for i in range(r - 1, -1, -1):
    c, pv = piv[i]
    s = M[i][-1]
    for j in range(c + 1, len(P)):
        s -= M[i][j] * sol.get(j, 0)
    if s % pv:
        bad.append((P[c], pv, s))
        say('   x_%d: %s / %s NOT divisible  (pivot %d bits, rhs %d bits)'
            % (P[c], 'rhs', 'pivot', abs(pv).bit_length(), abs(s).bit_length()))
        sol[c] = 0
    else:
        sol[c] = s // pv
        say('   x_%d = rhs/pivot OK (%d bits)' % (P[c], abs(sol[c]).bit_length()))
say('\nnumber of blocking divisibilities: %d' % len(bad))
for u, pv, s in bad:
    import math
    g = math.gcd(int(pv), int(s))
    say('   x_%d: gcd(pivot, rhs) = %d bits; pivot/gcd = %d bits' %
        (u, g.bit_length(), (abs(pv) // g).bit_length()))
say('DONE')
