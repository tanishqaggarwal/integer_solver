"""S11 step 11: the REFINED cost of every frame-2 drop set.

sacrifice2.py costed a drop set S as |union of equations touched|.  That is a bad
overcount: an equation containing several nonzero atoms fails only if their
COMBINATION is nonzero.  The delivered witness touches 12 equations with 7 atoms
and fails only 7.

For each S subset of the seven frame-2 checks: E = equations touched, M = the
|E| x |S| coefficient matrix.  The atoms in S may take any values A (the others
are 0), so an equation is satisfied iff M[r].A = 0.  Maximise the number of zero
rows over A != 0 -- i.e. find the largest subset of rows of rank < |S|.
"""
import os, sys, itertools, collections, json
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, ad
from fractions import Fraction
P = ad.P
SEVEN = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
eqsets = {a: set(L.atom2eq[a]) for a in SEVEN}
ALL = sorted(set().union(*eqsets.values()))
print(f'the seven residual atoms touch {len(ALL)} equations')

# coefficient of atom a in equation e
coef = {}
for e in ALL:
    _m, _sq, co = L.eq_atoms[e]
    coef[e] = {a: c for a, c in co.items() if a in SEVEN}

def rank_q(rows):
    """rank over Q of a list of integer vectors."""
    M = [list(map(Fraction, r)) for r in rows]
    n = len(M); m = len(M[0]) if n else 0
    r_ = 0
    for j in range(m):
        k = next((i for i in range(r_, n) if M[i][j] != 0), None)
        if k is None: continue
        M[r_], M[k] = M[k], M[r_]
        pv = M[r_][j]
        M[r_] = [x / pv for x in M[r_]]
        for i in range(n):
            if i != r_ and M[i][j] != 0:
                f = M[i][j]
                M[i] = [x - f * y for x, y in zip(M[i], M[r_])]
        r_ += 1
    return r_

best = None
results = []
for k in range(1, 8):
    for S in itertools.combinations(SEVEN, k):
        Sl = list(S)
        E = sorted(set().union(*[eqsets[a] for a in S]))
        rows = [[coef[e].get(a, 0) for a in Sl] for e in E]
        # equations whose row is all-zero are satisfied for free
        nzrows = [(e, r) for e, r in zip(E, rows) if any(r)]
        free_sat = len(E) - len(nzrows)
        # largest subset of nonzero rows with rank < k  (so a nonzero A kills them)
        bestsat = 0
        if k == 1:
            bestsat = 0
        else:
            idx = list(range(len(nzrows)))
            # greedy exact for small sizes: try every subset of rows up to rank k-1
            # by building maximal rank-(k-1) row subspaces from row combinations
            for combo in itertools.combinations(idx, min(k - 1, len(idx))):
                base = [nzrows[i][1] for i in combo]
                if not base: continue
                rb = rank_q(base)
                if rb >= k: continue
                cnt = 0
                for i in idx:
                    if rank_q(base + [nzrows[i][1]]) == rb: cnt += 1
                bestsat = max(bestsat, cnt)
        fail = len(E) - free_sat - bestsat
        results.append((fail, len(E), k, S))
        if best is None or fail < best[0]: best = (fail, len(E), k, S)
results.sort()
print(f'\n{"fail":>5} {"|E|":>4} {"|S|":>4}  drop set')
seen = set()
for fail, ne, k, S in results[:14]:
    print(f'{fail:>5} {ne:>4} {k:>4}  {S}')
print(f'\nBEST: drop {best[3]} -> {best[0]} failing equations '
      f'(score {L.NEQ - best[0]})')
print(f'delivered witness drops all seven: '
      f'{[r for r in results if r[3] == tuple(SEVEN)]}')
