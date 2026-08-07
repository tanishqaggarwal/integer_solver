"""WR step 12: EXACT max-satisfiable for the wire placement.

S = {a37694} plus compensators.  Region R(S) = union of their equations.  With
every atom of S treated as an INDEPENDENT free rational (a relaxation -- reality
has far fewer degrees of freedom), maximise #{e in R : row_e . v = 0} subject to
v[a37694] != 0.  Any failing count below this is impossible.
"""
import os, sys, collections, random, itertools, json
from fractions import Fraction
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
P = ad.P
ROOT = 37694


def region(S):
    R = set()
    for a in S:
        R |= set(L.atom2eq[a])
    return R


def maxsat(S, root=ROOT, tries=3000, seed=1, verbose=False):
    Sl = sorted(S)
    si = {a: i for i, a in enumerate(Sl)}
    n = len(Sl)
    Rl = sorted(region(S))
    rows = [[L.eq_atoms[e][2].get(a, 0) for a in Sl] for e in Rl]
    ri = si[root]
    rnd = random.Random(seed)
    best, bv, bset = -1, None, None
    for t in range(tries):
        order = list(range(len(Rl)))
        rnd.shuffle(order)
        M, piv = [], []
        for k in order:
            rr = [Fraction(x) for x in rows[k]]
            for i, pc in enumerate(piv):
                if rr[pc]:
                    f = rr[pc]
                    rr = [rr[j] - f * M[i][j] for j in range(n)]
            nzc = [j for j in range(n) if rr[j]]
            if not nzc:
                continue
            if nzc == [ri]:
                continue
            pc = next((j for j in nzc if j != ri), None)
            inv = rr[pc]
            rr = [x / inv for x in rr]
            M.append(rr); piv.append(pc)
        v = [Fraction(0)] * n
        v[ri] = Fraction(1)
        free = [j for j in range(n) if j not in piv and j != ri]
        for j in free:
            v[j] = Fraction(rnd.randrange(-5, 6))
        for i in reversed(range(len(piv))):
            pc = piv[i]
            s = sum(M[i][j] * v[j] for j in range(n) if j != pc and M[i][j])
            v[pc] = -s
        assert v[ri] == 1
        sat = [Rl[k] for k in range(len(Rl))
               if sum(Fraction(rows[k][j]) * v[j] for j in range(n) if rows[k][j]) == 0]
        if len(sat) > best:
            best, bv, bset = len(sat), v, sat
            if verbose:
                print(f'    try{t}: satisfied {best}/{len(Rl)} '
                      f'-> failing {len(Rl)-best}', flush=True)
    return best, len(Rl), bset, dict(zip(Sl, bv))


if __name__ == '__main__':
    print('=== S = {a37694} ===')
    b, tot, _, _ = maxsat({ROOT}, tries=5)
    print(f'   region {tot}, max satisfied {b}, FAILING {tot-b}')

    print('=== S = {a37694, a39417} ===')
    b, tot, _, _ = maxsat({ROOT, 39417}, tries=20)
    print(f'   region {tot}, max satisfied {b}, FAILING {tot-b}')

    # everything in the 13 equations
    S = set()
    for e in region({ROOT, 39417}):
        S |= set(L.eq_atoms[e][2])
    print(f'\n=== S = all atoms of those 13 equations ({len(S)} atoms) ===')
    b, tot, sat, _ = maxsat(S, tries=int(sys.argv[1]) if len(sys.argv) > 1 else 400,
                            verbose=True)
    print(f'   region {tot}, max satisfied {b}, FAILING {tot-b}')

    # one more closure round
    S2 = set(S)
    for e in region(S):
        S2 |= set(L.eq_atoms[e][2])
    print(f'\n=== closure round 2: S = {len(S2)} atoms ===')
    b2, tot2, _, _ = maxsat(S2, tries=int(sys.argv[2]) if len(sys.argv) > 2 else 120,
                            verbose=True)
    print(f'   region {tot2}, max satisfied {b2}, FAILING {tot2-b2}')
