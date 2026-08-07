"""W round 2, item 1+2: EXACT enumeration of every cocircuit of size <= 6 on K+ = 40.
NO window, NO degenerate skips -- this replaces round 1's heuristic s=3..6 search.

Method.  Fix an information set I (32 rows, identity block); C[j] is the dual basis.  Any
codeword c with |supp(c)| <= 6 vanishes on I\T for T = supp(c) cap I, |T| = s <= 6, so
c = sum_{j in T} y_j C[j] for some y in P^{s-1}.  Outside I there are NOUT columns and c may
be nonzero on at most 6-s of them.  Therefore:

    ** among ANY 7-s outside columns, at least one vanishes at c. **

That is an exact case split.  Branch on it: pick 7-s columns whose functional is not
identically zero on the current subspace V, recurse into V cap ker(phi_k) for each (each cut
drops dim by exactly 1).  Depth s-1, branching 7-s, so the tree is finite and COMPLETE.
If fewer than 7-s columns are nonzero on V, then every point of V has weight <= that count and
the whole subspace is recorded -- the case round 1 skipped as 'degenerate' is now handled.
"""
import sys, os, itertools, json, time
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup2 as S

P = (1 << 61) - 1
KN = S.KNOB; NK = len(KN); RANK = 32
SMAX = int(os.environ.get('WSMAX', '6'))
SMIN = int(os.environ.get('WSMIN', '1'))
ROWS = {e: [S.rows[e].get(u, 0) % P for u in KN] for e in S.SAT}

def gauss(M, n):
    M = [r[:] for r in M]; piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c]), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        inv = pow(M[r][c], P-2, P); M[r] = [x*inv % P for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(a-f*b) % P for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
        if r == len(M): break
    return M[:r], piv

I, basis = [], []
for e in S.SAT:
    if len(gauss(basis+[ROWS[e]], NK)[0]) > len(basis): basis = basis+[ROWS[e]]; I.append(e)
    if len(I) == RANK: break
assert len(I) == RANK
OUT = [e for e in S.SAT if e not in I]; NOUT = len(OUT)
AI = [ROWS[e] for e in I]
aug = [AI[i] + [1 if k == i else 0 for k in range(RANK)] for i in range(RANK)]
R, piv = gauss(aug, NK)
C = [[0]*NK for _ in range(RANK)]
for i, c in enumerate(piv):
    for j in range(RANK): C[j][c] = R[i][NK+j]
for i, e in enumerate(I):
    for j in range(RANK):
        assert sum(a*b for a, b in zip(ROWS[e], C[j])) % P == (1 if i == j else 0)
print('information set %d rows, outside %d, dual basis verified' % (RANK, NOUT), flush=True)
# GT[j][k]
G = [[sum(a*b for a, b in zip(ROWS[e], C[j])) % P for e in OUT] for j in range(RANK)]
GC = [[G[j][k] for j in range(RANK)] for k in range(NOUT)]   # column-major

found = set(); subspace_hits = 0; t0 = time.time(); leaves = 0

def solve_T(Tidx):
    global subspace_hits, leaves
    s = len(Tidx); budget = 6 - s; br = budget + 1        # = 7-s
    col = [[GC[k][j] for j in Tidx] for k in range(NOUT)]  # NOUT x s
    def rec(V):
        global subspace_hits, leaves
        d = len(V)
        if d == 1:
            leaves += 1
            y = V[0]; nz = []
            for k in range(NOUT):
                ck = col[k]; v = 0
                for a in range(s): v += y[a]*ck[a]
                if v % P:
                    nz.append(k)
                    if len(nz) > budget: return
            found.add(frozenset([I[j] for j in Tidx] + [OUT[k] for k in nz])); return
        live = []
        for k in range(NOUT):
            ck = col[k]
            for b in V:
                v = 0
                for a in range(s): v += b[a]*ck[a]
                if v % P: live.append(k); break
            if len(live) == br: break
        if len(live) < br:
            subspace_hits += 1
            found.add(frozenset([I[j] for j in Tidx] + [OUT[k] for k in live])); return
        for k in live:
            ck = col[k]
            vals = [sum(b[a]*ck[a] for a in range(s)) % P for b in V]
            i0 = next(i for i, v in enumerate(vals) if v)
            inv = pow(vals[i0], P-2, P)
            V2 = [[(V[i][a] - vals[i]*inv % P * V[i0][a]) % P for a in range(s)]
                  for i in range(d) if i != i0]
            rec(V2)
    rec([[1 if a == b else 0 for a in range(s)] for b in range(s)])

for s in range(SMIN, SMAX+1):
    n0 = len(found); ts = time.time()
    for Tidx in itertools.combinations(range(RANK), s): solve_T(Tidx)
    print('  s=%d EXACT: %d new supports (total %d)  %.0fs  [leaves %d, subspace-nodes %d]'
          % (s, len(found)-n0, len(found), time.time()-ts, leaves, subspace_hits), flush=True)
print('\nEXACT enumeration complete for s = %d..%d.  supports: %d  elapsed %.0fs'
      % (SMIN, SMAX, len(found), time.time()-t0))
print('positive-dimensional (whole-subspace) nodes recorded rather than skipped: %d' % subspace_hits)
json.dump({'supports': [sorted(map(str, T)) for T in sorted(found, key=lambda t: (len(t), sorted(map(str, t))))],
           'exact': True, 'smin': SMIN, 'smax': SMAX, 'subspace_hits': subspace_hits},
          open('w_cocirc3_raw_s%d_%d.json' % (SMIN, SMAX), 'w'), indent=1)
