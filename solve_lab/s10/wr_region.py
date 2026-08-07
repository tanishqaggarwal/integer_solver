"""WR step 7: EXACT price of the placement S = {a37694, ...}.

With a37694 != 0 forced, the failing count over the region R(S) = U eqs(a) obeys
    #satisfied <= (max rows of R lying in one hyperplane not containing e_37694)
              <= |S| - 1               (rows in general position)
so   #failing >= |R(S)| - |S| + 1.
Grow S greedily by the atom that drags in the fewest new equations and watch the
bound.  Then compute the ACTUAL max-satisfiable by exact linear algebra.
"""
import os, sys, collections, random, itertools
from fractions import Fraction
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
ROOTATOM = 37694


def region(S):
    R = set()
    for a in S:
        R |= set(L.atom2eq[a])
    return R


def grow(S0, steps=40, verbose=True):
    S = set(S0)
    R = region(S)
    hist = []
    for it in range(steps):
        bound = len(R) - len(S) + 1
        hist.append((len(S), len(R), bound))
        if verbose:
            print(f'  |S|={len(S):<4} |R|={len(R):<5} bound_failing>={bound}')
        # candidate atoms: any atom appearing in R, not already in S
        cand = set()
        for e in R:
            cand |= set(L.eq_atoms[e][2])
        cand -= S
        best = None
        for a in cand:
            new = len(set(L.atom2eq[a]) - R)
            if best is None or new < best[0]:
                best = (new, a)
        if best is None:
            break
        new, a = best
        S.add(a)
        R |= set(L.atom2eq[a])
        if verbose and it < 30:
            print(f'      + a{a} (drags {new} new eqs, total eqs {len(L.atom2eq[a])})')
    return S, R, hist


def maxsat_exact(S, R, root=ROOTATOM, tries=4000, seed=1):
    """Randomised hyperplane search: max #rows of R vanishing at some v with
    v[root] = 1.  Returns (best_count, best_v)."""
    Sl = sorted(S)
    si = {a: i for i, a in enumerate(Sl)}
    n = len(Sl)
    Rl = sorted(R)
    rows = []
    for e in Rl:
        co = L.eq_atoms[e][2]
        rows.append([co.get(a, 0) for a in Sl])
    ri = si[root]
    rnd = random.Random(seed)
    best, bv = -1, None
    for t in range(tries):
        order = list(range(len(Rl)))
        rnd.shuffle(order)
        # build constraint system incrementally: keep rows that don't force v[root]=0
        M = []            # reduced rows (list of Fraction lists) with pivots
        piv = []
        chosen = []
        for k in order:
            r = [Fraction(x) for x in rows[k]]
            rr = r[:]
            for i, pc in enumerate(piv):
                if rr[pc]:
                    f = rr[pc]
                    rr = [rr[j] - f * M[i][j] for j in range(n)]
            nzc = [j for j in range(n) if rr[j]]
            if not nzc:
                chosen.append(k)
                continue
            if nzc == [ri]:
                continue            # would force v[root] = 0
            # normalise on a pivot that is not the root column when possible
            pc = next((j for j in nzc if j != ri), nzc[0])
            inv = rr[pc]
            rr = [x / inv for x in rr]
            M.append(rr); piv.append(pc); chosen.append(k)
        # solve: free vars -> set v[root] = 1, others 0 where free
        v = [Fraction(0)] * n
        v[ri] = Fraction(1)
        for i in reversed(range(len(piv))):
            pc = piv[i]
            s = sum(M[i][j] * v[j] for j in range(n) if j != pc)
            v[pc] = -s
        if v[ri] != 1:
            continue
        cnt = 0
        for k in range(len(Rl)):
            if sum(Fraction(rows[k][j]) * v[j] for j in range(n) if rows[k][j]) == 0:
                cnt += 1
        if cnt > best:
            best, bv = cnt, v
    return best, bv


if __name__ == '__main__':
    print('=== placement S = {a37694} ===')
    S0 = {ROOTATOM}
    print(f'|R| = {len(region(S0))}, bound_failing >= {len(region(S0)) - 1 + 1}')
    S, R, hist = grow(S0, steps=int(sys.argv[1]) if len(sys.argv) > 1 else 25)
    print(f'\nfinal |S|={len(S)} |R|={len(R)} bound={len(R)-len(S)+1}')
    print(f'best bound over the growth: {min(h[2] for h in hist)}')

    print('\n=== placement S = {a37694, a39417} (uniform wire shift) ===')
    S1 = {ROOTATOM, 39417}
    R1 = region(S1)
    print(f'|R| = {len(R1)}, bound_failing >= {len(R1) - len(S1) + 1}')
    S, R, hist = grow(S1, steps=int(sys.argv[1]) if len(sys.argv) > 1 else 25)
    print(f'final |S|={len(S)} |R|={len(R)} bound={len(R)-len(S)+1}')
    print(f'best bound over the growth: {min(h[2] for h in hist)}')
