"""Exact rational structure of O's 175x34 frame-B system.

For each P (subset of the 7 failing rows we insist on buying):
   A_P   = affine solution set of the rows in P
   Z_P   = SAT rows whose affine form vanishes IDENTICALLY on A_P   (free -- stay satisfied)
   Q_P   = SAT rows that do not                                     (each vanishes on a proper subset)
   dim A_P
If dim A_P == 0 the answer is immediate: A_P is a single rational point and minbreak(P) is just
the number of SAT rows nonzero there (or +inf if that point is not an integer point).
"""
import sys, os, itertools, json
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S

KN = S.KNOB
NK = len(KN)
kidx = {u: i for i, u in enumerate(KN)}


def mat(names):
    return [[Fraction(S.rows[e].get(u, 0)) for u in KN] + [Fraction(S.rhs[e])] for e in names]


def rref(M):
    M = [r[:] for r in M]
    m = len(M); n = len(M[0])
    piv = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        d = M[r][c]
        M[r] = [x / d for x in M[r]]
        for i in range(m):
            if i != r and M[i][c] != 0:
                f = M[i][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == m: break
    return M[:r], piv


print('coef bits: max %d' % max(abs(v).bit_length() for e in S.names for v in S.rows[e].values() if v))
print('rhs  bits: max %d' % max(abs(S.rhs[e]).bit_length() for e in S.names if S.rhs[e]))

full, fp = rref(mat(S.names))
print('rank of ALL 175 rows (incl. rhs col) = %d ; consistent-cols rank=%d'
      % (len(full), len([c for c in fp if c < NK])))
sat, sp = rref(mat(S.SAT))
print('rank of the 168+S satisfied rows: %d, rhs is pivot: %s' % (len(sat), NK in sp))
print('dim of solution set of SAT alone = %d' % (NK - len([c for c in sp if c < NK])))

out = {}
for k in range(1, 8):
    for P in itertools.combinations(S.FAIL, k):
        R, pv = rref(mat(list(P)))
        if NK in pv:
            out[str(list(P))] = {'A_P': 'EMPTY over Q'}
            print('P=%-34s A_P EMPTY over Q' % (list(P),)); continue
        dim = NK - len(pv)
        # Z_P: rows whose affine form is in the row space of P's system (identically zero on A_P)
        Z = []
        base = [r[:] for r in R]
        for e in S.SAT:
            row = [Fraction(S.rows[e].get(u, 0)) for u in KN] + [Fraction(S.rhs[e])]
            v = row[:]
            for i, c in enumerate(pv):
                if v[c] != 0:
                    f = v[c]; v = [a - f * b for a, b in zip(v, base[i])]
            if all(x == 0 for x in v): Z.append(e)
        Q = [e for e in S.SAT if e not in Z]
        out[str(list(P))] = {'dimA': dim, 'nZ': len(Z), 'nQ': len(Q), 'Q': [str(x) for x in Q]}
        print('P=%-34s dim A_P=%2d  free-SAT |Z|=%3d  contested |Q|=%2d  Q=%s'
              % (list(P), dim, len(Z), len(Q), [str(x) for x in Q]))
json.dump(out, open('w_struct.json', 'w'), indent=1)
