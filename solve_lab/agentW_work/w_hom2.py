"""W: preconditions for the cocircuit argument on K+ = 40 -- verify BEFORE porting."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup2 as S
from fractions import Fraction
nzr = [e for e in S.SAT if S.rhs[e] != 0]
print('SAT rows %d ; with NONZERO rhs: %d  %s' % (len(S.SAT), len(nzr), nzr[:6]))
print('FAIL rows %d ; rhs: %s' % (len(S.FAIL), {e: S.rhs[e] != 0 for e in S.FAIL}))
def rank_of(rows, n):
    basis = []
    for v in rows:
        w = v[:]
        for bp, bv in basis:
            if w[bp] != 0:
                f = w[bp]; w = [a - f*b for a, b in zip(w, bv)]
        p = next((i for i, x in enumerate(w) if x != 0), None)
        if p is not None:
            dd = w[p]; basis.append((p, [x/dd for x in w]))
    return len(basis)
KN = S.KNOB
V = lambda e: [Fraction(S.rows[e].get(u, 0)) for u in KN]
r = rank_of([V(e) for e in S.SAT], len(KN))
ra = rank_of([V(e) + [Fraction(S.rhs[e])] for e in S.SAT], len(KN)+1)
print('rank(A_SAT) = %d ; rank([A_SAT|b]) = %d  -> homogeneous: %s' % (r, ra, r == ra))
rall = rank_of([V(e) for e in S.SAT+S.FAIL], len(KN))
rallb = rank_of([V(e)+[Fraction(S.rhs[e])] for e in S.SAT+S.FAIL], len(KN)+1)
print('rank(A_ALL) = %d ; rank([A_ALL|b]) = %d  -> Q-consistent: %s' % (rall, rallb, rall == rallb))
