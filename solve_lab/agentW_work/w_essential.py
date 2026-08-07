"""Which SAT rows can breaking possibly help?

SAT rows are HOMOGENEOUS (rhs 0), so the set of admissible knob deltas keeping a kept-set
KEEP satisfied is exactly ker_Z(A_KEEP) = Z^34 cap ker_Q(A_KEEP).  Hence if
rank_Q(A_KEEP) == rank_Q(A_SAT), the integer feasible set is IDENTICAL -- breaking those rows
buys literally nothing, at any budget.  So only rank-reducing deletions matter.

A row r is ESSENTIAL (its single deletion drops the rank) iff e_r lies in colspace(A_SAT),
i.e. iff the rref of A_SAT^T has the unit vector e_r as a row.
"""
import sys, os, json
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S

KN = S.KNOB
print('SAT rows with nonzero rhs:', [e for e in S.SAT if S.rhs[e] != 0])
print('FAIL rows rhs nonzero:', all(S.rhs[e] != 0 for e in S.FAIL))

def rref(M):
    M = [r[:] for r in M]; m = len(M); n = len(M[0]); piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        d = M[r][c]; M[r] = [x / d for x in M[r]]
        for i in range(m):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == m: break
    return M[:r], piv

# A_SAT^T : 34 x 169 (rows = knob columns of A_SAT)
AT = [[Fraction(S.rows[e].get(u, 0)) for e in S.SAT] for u in KN]
R, piv = rref(AT)
print('rank(A_SAT) = %d  (rref of A_SAT^T has %d pivot rows)' % (len(R), len(R)))
ess = []
for i, c in enumerate(piv):
    row = R[i]
    if all((x == 0) for j, x in enumerate(row) if j != c) and row[c] == 1:
        ess.append(S.SAT[c])
print('ESSENTIAL SAT rows (deleting one drops rank 26 -> 25): %d  -> %s' % (len(ess), ess))
red = [e for e in S.SAT if e not in ess]
print('REDUNDANT SAT rows: %d  -- deleting ANY subset of these alone cannot change the integer feasible set' % len(red))
json.dump({'essential': [str(x) for x in ess], 'nred': len(red)}, open('w_essential.json','w'), indent=1)
