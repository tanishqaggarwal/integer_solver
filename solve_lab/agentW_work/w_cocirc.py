"""Enumerate every low-cost direction: all SAT-row sets of size <= 6 that a knob delta can break.

Why this is the whole remaining question.  SAT rows are homogeneous, so the admissible deltas
keeping KEEP satisfied are ker_Z(A_KEEP) = Z^34 cap ker_Q(A_KEEP).  Hence if breaking B does
not drop rank_Q, the feasible set is unchanged and breaking B buys nothing.  So the only break
sets that can ever help are those containing a COCIRCUIT (a minimal rank-dropping set).  Every
assignment with cost <= 6 breaks such a set.  Enumerate all of them with size <= 6 and the
frame-B budget is closed at every j <= 7, not just j <= 3.

Method: the map c -> (row_e . c)_e is a rank-26 code of length 169.  Fix an information set I
(26 rows, identity block).  A codeword of weight <= 6 has <= 6 nonzeros on I, so its support
on I is some T, |T| = s <= 6, and it kills all but <= 6-s of the other 143 columns.  Among any
5 of those columns at least s-1 are killed, and s-1 killed columns pin the codeword up to
scale.  So enumerating T and the C(5, s-1) subsets of a 5-column window is complete, except
where every such subset is degenerate -- checked for, counted, and reported.

Arithmetic is mod a 61-bit prime for the search; every hit is re-verified exactly over Q.
"""
import sys, os, itertools, json, time
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup as S

P = (1 << 61) - 1
KN = S.KNOB
NK = len(KN)
ROWS = {e: [S.rows[e].get(u, 0) % P for u in KN] for e in S.SAT}


def gauss(M, n):
    """rref of M (list of rows length n) mod P; returns (rref rows, pivots)."""
    M = [r[:] for r in M]; piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        inv = pow(M[r][c], P - 2, P)
        M[r] = [x * inv % P for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(a - f * b) % P for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    return M[:r], piv


# information set
I, basis = [], []
for e in S.SAT:
    trial = basis + [ROWS[e]]
    if len(gauss(trial, NK)[0]) > len(basis):
        basis = trial; I.append(e)
    if len(I) == 26: break
assert len(I) == 26, len(I)
OUT = [e for e in S.SAT if e not in I]
print('information set: 26 rows, %d outside' % len(OUT), flush=True)

# dual vectors: c_j with row_i . c_j = delta_ij for i in I
AI = [ROWS[e] for e in I]
aug = [AI[i] + [1 if k == i else 0 for k in range(26)] for i in range(26)]
R, piv = gauss(aug, NK)
C = [[0] * NK for _ in range(26)]      # C[j] = c_j
for i, c in enumerate(piv):
    for j in range(26):
        C[j][c] = R[i][NK + j]
for i, e in enumerate(I):
    for j in range(26):
        v = sum(a * b for a, b in zip(ROWS[e], C[j])) % P
        assert v == (1 if i == j else 0), (i, j, v)
print('dual basis verified (identity block)', flush=True)

G = [[sum(a * b for a, b in zip(ROWS[e], C[j])) % P for e in OUT] for j in range(26)]

NOUT = len(OUT)
found = set(); degenerate = 0; t0 = time.time()

def record(Tidx, y):
    nz = [k for k in range(NOUT) if sum(y[a] * G[j][k] for a, j in enumerate(Tidx)) % P]
    if len(nz) <= 6 - len(Tidx):
        found.add(frozenset([I[j] for j in Tidx] + [OUT[k] for k in nz])); return 1
    return 0

# s = 1 and s = 2 are handled EXACTLY (no window, no degeneracy):
#   s=1: the codeword is the basis codeword itself.
#   s=2: y = (1, lam); a column k is killed iff lam = -G[i,k]/G[j,k], so the best lam is the
#        most frequent such ratio -- computed over ALL 142 columns, so this is exhaustive.
cnt = 0
for j in range(26): cnt += record((j,), [1])
print('  s=1 (exact): %d hits, %d supports, %.0fs' % (cnt, len(found), time.time()-t0), flush=True)
cnt = 0
for i, j in itertools.combinations(range(26), 2):
    from collections import Counter
    ratios = Counter(); always = 0
    for k in range(NOUT):
        a, b = G[i][k], G[j][k]
        if a == 0 and b == 0: always += 1
        elif b: ratios[(-a) * pow(b, P-2, P) % P] += 1
    for lam, m in ratios.most_common(4):
        if always + m >= NOUT - 4: cnt += record((i, j), [1, lam])
print('  s=2 (exact): %d hits, %d supports, %.0fs' % (cnt, len(found), time.time()-t0), flush=True)

# s >= 3: window over the first WIN columns on which T is nonzero (adaptive, avoids the
# all-zero-column degeneracy).  Complete only if some (s-1)-subset of the window is both
# killed and independent; degenerate subsets are counted and reported.
for s in range(3, 7):
    WIN = 8 if s <= 4 else 6
    cnt = 0
    for Tidx in itertools.combinations(range(26), s):
        win = []
        for k in range(NOUT):
            if any(G[j][k] for j in Tidx):
                win.append(k)
                if len(win) == WIN: break
        for Z in itertools.combinations(win, s - 1):
            M = [[G[j][k] for j in Tidx] for k in Z]
            Rr, pv = gauss(M, s)
            free = [c for c in range(s) if c not in pv]
            if len(free) != 1:
                degenerate += 1; continue
            f = free[0]
            y = [0]*s; y[f] = 1
            for i2, c in enumerate(pv): y[c] = (-Rr[i2][f]) % P
            cnt += record(Tidx, y)
    print('  s=%d (window %d): %d hits, %d supports, %.0fs' % (s, WIN, cnt, len(found), time.time()-t0), flush=True)
print('degenerate (skipped) subsets: %d' % degenerate)

# exact verification over Q of every candidate
def rankQ(rows):
    b = []
    for e in rows:
        w = [Fraction(S.rows[e].get(u, 0)) for u in KN]
        for bp, bv in b:
            if w[bp] != 0:
                f = w[bp]; w = [a - f * x for a, x in zip(w, bv)]
        p = next((i for i, x in enumerate(w) if x != 0), None)
        if p is not None:
            d = w[p]; b.append((p, [x / d for x in w]))
    return len(b)


ver = []
for T in sorted(found, key=lambda t: (len(t), sorted(map(str, t)))):
    keep = [e for e in S.SAT if e not in T]
    if rankQ(keep) < 26:
        ver.append(sorted(T, key=str))
print('\nCOCIRCUIT-CONTAINING SETS of size <= 6, exactly verified over Q: %d' % len(ver))
for T in ver:
    print('   size %d : %s' % (len(T), [str(x) for x in T]))
json.dump([[str(x) for x in T] for T in ver], open('w_cocirc.json', 'w'), indent=1)
print('elapsed %.0fs' % (time.time() - t0))
