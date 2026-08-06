"""WR step 15: uniform w=1 PLUS three copy-atom knobs.

d = c*1 + t1*e_13720 + t2*e_12752 + t3*e_18306 makes
   a37691 = t1, a37692 = t2, a37693 = t3, a37694 = c
nonzero and everything else in the identity system unchanged (up to the rows that
carry those three coordinates).  Enumerate which rows can be zeroed and price the
best choice exactly, over ALL triples of wire coordinates.
"""
import os, sys, collections, itertools, json
from fractions import Fraction
HERE = '/home/user/integer_solver/solve_lab/s10'
LAB = '/home/user/integer_solver/solve_lab'
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
sys.path.insert(0, HERE)
import lib as L, tools as T, ad
import wr_rows as R
P = ad.P
WIRE, widx, rows, RE = R.WIRE, R.widx, R.rows, R.RE
N = len(WIRE)
rowsum = {e: sum(rows[e].values()) for e in RE}
BAD12 = [e for e in RE if rowsum[e]]
print(f'rows with nonzero row-sum (broken by a uniform shift): {len(BAD12)} -> {BAD12}')

# choose k coordinates and solve to zero k of the BAD12 rows; count total broken.
def price(coords, target, c=1):
    """coords: list of wire indices; target: list of rows to zero (same length).
    Solve A t = -c*rowsum(target); return (#broken rows, t)."""
    k = len(coords)
    A = [[Fraction(rows[e].get(j, 0)) for j in coords] for e in target]
    b = [Fraction(-c * rowsum[e]) for e in target]
    M = [A[i] + [b[i]] for i in range(k)]
    piv = []
    r = 0
    for col in range(k):
        s = next((i for i in range(r, k) if M[i][col]), None)
        if s is None:
            continue
        M[r], M[s] = M[s], M[r]
        inv = M[r][col]
        M[r] = [x / inv for x in M[r]]
        for i in range(k):
            if i != r and M[i][col]:
                f = M[i][col]
                M[i] = [M[i][j] - f * M[r][j] for j in range(k + 1)]
        piv.append(col); r += 1
    for i in range(r, k):
        if M[i][k]:
            return None                      # inconsistent
    t = [Fraction(0)] * k
    for i, col in enumerate(piv):
        t[col] = M[i][k]
    broken = []
    for e in RE:
        s = Fraction(c * rowsum[e])
        for i, j in enumerate(coords):
            cc = rows[e].get(j, 0)
            if cc and t[i]:
                s += cc * t[i]
        if s:
            broken.append(e)
    return broken, t


if __name__ == '__main__':
    TRI = [widx[u] for u in (13720, 12752, 18306)]
    best = None
    for k in (1, 2, 3):
        loc = None
        for tgt in itertools.combinations(BAD12, k):
            for coords in itertools.combinations(TRI, k):
                r = price(list(coords), list(tgt))
                if r is None:
                    continue
                broken, t = r
                if loc is None or len(broken) < len(loc[0]):
                    loc = (broken, t, coords, tgt)
        if loc:
            print(f'k={k} knobs from (13720,12752,18306): best {len(loc[0])} broken rows '
                  f'-> {loc[0]}')
            if best is None or len(loc[0]) < len(best[0]):
                best = loc

    print('\nnow ALL coordinate triples (220 choose 3 with pruning):')
    # k = 1: every coordinate, every target row
    b1 = None
    for j in range(N):
        for e in BAD12:
            r = price([j], [e])
            if r is None:
                continue
            broken, t = r
            if b1 is None or len(broken) < len(b1[0]):
                b1 = (broken, t, (j,), (e,))
    print(f'  k=1 best: {len(b1[0])} broken rows, coord x_{WIRE[b1[2][0]]}, row {b1[3][0]}')
    print(f'      broken: {b1[0]}')
    # k = 2 restricted to the 40 coordinates that appear in the BAD12 rows
    cand = sorted(set().union(*[set(rows[e]) for e in BAD12]))
    print(f'  coordinates occurring in the 12 broken rows: {len(cand)}')
    b2 = None
    for coords in itertools.combinations(cand, 2):
        for tgt in itertools.combinations(BAD12, 2):
            r = price(list(coords), list(tgt))
            if r is None:
                continue
            broken, t = r
            if b2 is None or len(broken) < len(b2[0]):
                b2 = (broken, t, coords, tgt)
    print(f'  k=2 best: {len(b2[0])} broken rows, coords '
          f'{[WIRE[j] for j in b2[2]]}, rows {b2[3]}')
    print(f'      broken: {b2[0]}')
    b3 = None
    for coords in itertools.combinations(cand, 3):
        for tgt in itertools.combinations(BAD12, 3):
            r = price(list(coords), list(tgt))
            if r is None:
                continue
            broken, t = r
            if b3 is None or len(broken) < len(b3[0]):
                b3 = (broken, t, coords, tgt)
                if len(broken) <= 5:
                    break
    print(f'  k=3 best: {len(b3[0])} broken rows, coords '
          f'{[WIRE[j] for j in b3[2]]}, rows {b3[3]}')
    print(f'      broken: {b3[0]}')
