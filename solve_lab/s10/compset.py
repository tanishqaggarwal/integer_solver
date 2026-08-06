"""S11 step 24: cheapest COMPENSATOR SET that lets all 12 equations hold with A != 0.

We need  M7*A + Mk*B = 0  with A admissible and B free, i.e. M7*A in colspace(Mk).
Search subsets of the cheap compensators; cost = |union of their equations OUTSIDE
the twelve|.  Anything under 7 beats the deliverable.
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
ESET = set(E)
col = {}
for e in E:
    m, sq, co = L.eq_atoms[e]
    for a, c in co.items(): col.setdefault(a, {})[e] = c
def vec(a): return [col.get(a, {}).get(e, 0) for e in E]
M7 = [vec(a) for a in SEVEN]                     # 7 columns, each length 12

cands = []
for a in col:
    if a in SEVEN: continue
    out = sorted(set(L.atom2eq[a]) - ESET)
    cands.append((len(out), a, out))
cands.sort()
print(f'{len(cands)} compensator candidates; cheapest: '
      f'{[(a, n) for n, a, _ in cands[:10]]}')

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

r7 = rank_q(M7)
print(f'rank(M7) = {r7} of 7 columns, 12 rows')
pool = [(n, a, out) for n, a, out in cands if n <= 8][:18]
print(f'searching subsets of the {len(pool)} cheapest (outside <= 8)')
best = None
for k in range(1, 8):
    found = False
    for S in itertools.combinations(pool, k):
        eqs = set()
        for n, a, out in S: eqs |= set(out)
        if best and len(eqs) >= best[0]: continue
        cols_ = M7 + [vec(a) for n, a, out in S]
        if rank_q(cols_) < 7 + k:                 # a kernel exists
            # does the kernel have a nonzero A-part?
            found = True
            if best is None or len(eqs) < best[0]:
                best = (len(eqs), [a for n, a, out in S], k)
                print(f'  k={k}: atoms {[a for n,a,_ in S]} -> kernel; '
                      f'cost {len(eqs)} equations outside', flush=True)
    if best and best[2] == k:
        pass
print(f'\nBEST compensator set: {best}')
if best:
    print(f'  -> all 12 equations could hold, at a cost of {best[0]} equations '
          f'elsewhere (score {L.NEQ - best[0]})')
