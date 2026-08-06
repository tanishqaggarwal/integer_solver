"""S10 step 12: EXACT solution of the residual optimisation.

Achievable set (derived, then verified numerically below).  With
  D  = x_7068 - x_2099           (both untouched by these knobs)
  A  = (a22229,a22230,a35758,a35759,a35760,a35761,a35762)
the knobs x_642,x_17325,x_9413,x_28730,x_29854,x_1329,x_31864,x_10903,x_9118,x_8731
realise EXACTLY

     A1 + 7376877*A7  ==  D   (mod 7376877*p)        <-- the only constraint
     A2..A6 free

Proof: a22229 = D - 7376877*t1, a35762 = t1 - p*t2  =>  A1 + 7376877*A7 = D - 7376877*p*t2.
       a22230 = t8 - p*t3 (free); a35758 = t4 - p*t5, a35759 = 5113045*t9 - t4 (free pair,
       gcd(5113045,p)=1); a35760 = t6 - p*t7, a35761 = t10 + t6 (free pair).

So: enumerate subsets S of the 12 touched equations, keep those for which
{A : M_S A = 0} contains a point of the achievable set.  Score = 39033 - 12 + |S|.
"""
import os, sys, json, itertools
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
LAB = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(LAB, 's9', 'eff'))
import lib as L

P = 2**256 - 2**32 - 977
MOD = 7376877 * P
BEST = os.path.join(LAB, 'best', 'new_instance_partial_39026.json')

v = L.load(BEST)
av = L.all_atom_values(v)
NZ = [22229, 22230, 35758, 35759, 35760, 35761, 35762]
EQS = sorted(L.eqs_of_atoms(NZ))
M = [[L.eq_atoms[i][2].get(a, 0) for a in NZ] for i in EQS]
D = v[7068] - v[2099]
A0 = [av[a] for a in NZ]
print('D =', D)
print('sanity: A1 + 7376877*A7 - D  mod 7376877p =',
      (A0[0] + 7376877 * A0[6] - D) % MOD)
sat0 = [i for i in EQS if L.eq_value(i, av) == 0]
print(f'base satisfied {len(sat0)}/{len(EQS)}: {sat0}  -> score {L.NEQ-len(EQS)+len(sat0)}')


# ---------- integer kernel via column-style HNF with unimodular tracking -----
def int_kernel(mat):
    """Z-basis (list of column vectors) of {x in Z^n : mat x = 0}."""
    m = len(mat); n = len(mat[0])
    A = [row[:] for row in mat]
    U = [[1 if i == j else 0 for j in range(n)] for i in range(n)]  # U[i][j]: col j
    piv_row = 0; piv_cols = []
    for r in range(m):
        # find columns (not yet pivot) with nonzero entry in row r
        while True:
            nz = [c for c in range(n) if c not in piv_cols and A[r][c] != 0]
            if len(nz) <= 1:
                break
            # euclidean reduce: use column with smallest |value| as pivot
            nz.sort(key=lambda c: abs(A[r][c]))
            p0 = nz[0]
            for c in nz[1:]:
                q = A[r][c] // A[r][p0]
                if q:
                    for i in range(m): A[i][c] -= q * A[i][p0]
                    for i in range(n): U[i][c] -= q * U[i][p0]
        nz = [c for c in range(n) if c not in piv_cols and A[r][c] != 0]
        if nz:
            piv_cols.append(nz[0])
    ker = [c for c in range(n) if c not in piv_cols]
    return [[U[i][c] for i in range(n)] for c in ker]


def solvable(S):
    """Is there an achievable A with M_S A = 0?  Return a witness A or None."""
    if not S:
        rows = []
    else:
        rows = [M[EQS.index(i)] for i in S]
    if rows:
        B = int_kernel(rows)          # basis vectors of ker_Z(M_S)
    else:
        B = [[1 if i == j else 0 for i in range(7)] for j in range(7)]
    if not B:
        return None
    # need y in Z^k with  sum_j y_j * (B_j[0] + 7376877*B_j[6]) == D  (mod MOD)
    c = [b[0] + 7376877 * b[6] for b in B]
    g = 0
    for x in c: g = __import__('math').gcd(g, x)
    g = __import__('math').gcd(g, MOD)
    if D % g:
        return None
    # construct a witness with extended gcd over the c_j
    import math
    # solve sum y_j c_j = D mod MOD  -> find combination giving gcd, then scale
    cur, coef = 0, [0] * len(c)
    for j, x in enumerate(c):
        if x == 0: continue
        if cur == 0:
            cur, coef = x, [0] * len(c); coef[j] = 1
        else:
            # extended gcd of cur and x
            a, b = cur, x
            old_r, r = a, b; old_s, s = 1, 0; old_t, t = 0, 1
            while r:
                q = old_r // r
                old_r, r = r, old_r - q * r
                old_s, s = s, old_s - q * s
                old_t, t = t, old_t - q * t
            coef = [old_s * z for z in coef]; coef[j] += old_t
            cur = old_r
    gg = math.gcd(cur if cur else 0, MOD)
    if gg == 0 or D % gg:
        return None
    mult = (D // gg) * pow((cur // gg) % (MOD // gg), -1, MOD // gg) % (MOD // gg)
    y = [mult * z for z in coef]
    A = [sum(y[j] * B[j][i] for j in range(len(B))) for i in range(7)]
    # verify
    for i in S:
        if sum(M[EQS.index(i)][k] * A[k] for k in range(7)) != 0:
            return None
    if (A[0] + 7376877 * A[6] - D) % MOD:
        return None
    return A


best = (len(sat0), tuple(sat0), None)
for size in range(12, 0, -1):
    found = []
    for S in itertools.combinations(EQS, size):
        A = solvable(list(S))
        if A is not None:
            found.append((S, A))
    print(f'  subsets of size {size:>2}: {len(found)} integer-solvable')
    if found and size > best[0]:
        best = (size, found[0][0], found[0][1])
        print(f'    -> IMPROVEMENT: satisfy {size} of 12  => score '
              f'{L.NEQ-len(EQS)+size}')
        print(f'       S = {found[0][0]}')
        print(f'       A = {[str(x)[:40] for x in found[0][1]]}')
        json.dump({'S': list(found[0][0]), 'A': [str(x) for x in found[0][1]],
                   'EQS': EQS, 'NZ': NZ, 'D': str(D)},
                  open(os.path.join(HERE, 'lattice_best.json'), 'w'))
        break
print(f'\nBEST: satisfy {best[0]} of the 12 -> score {L.NEQ-len(EQS)+best[0]}')
