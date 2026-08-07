"""W round 2, item 1: redo the cocircuit / general-break closure on the CORRECTED knob set
K+ = 40.  Preconditions re-verified in w_hom2.py: the 198 SAT rows are homogeneous and
rank(A_SAT) = 32, so admissible deltas for a kept-set are ker_Z(A_KEEP) and only RANK-DROPPING
deletions can ever help.  Every break-set of size <= 6 is a union of minimal cocircuits.

Change vs round 1 (besides K -> K+): the s=3..6 search is run from SEVERAL RANDOM INFORMATION
SETS and the supports unioned, so the marginal yield of each extra information set is itself a
measurement of how saturated the enumeration is.
"""
import sys, os, itertools, json, time, random
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import w_setup2 as S

P = (1 << 61) - 1
KN = S.KNOB; NK = len(KN)
ROWS = {e: [S.rows[e].get(u, 0) % P for u in KN] for e in S.SAT}
RANK = 32
NSETS = int(sys.argv[1]) if len(sys.argv) > 1 else 3

def gauss(M, n):
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

found = set(); degen_tot = 0; t0 = time.time(); yields = []
for trial in range(NSETS):
    order = list(S.SAT)
    if trial: random.Random(1000 + trial).shuffle(order)
    I, basis = [], []
    for e in order:
        if len(gauss(basis + [ROWS[e]], NK)[0]) > len(basis):
            basis = basis + [ROWS[e]]; I.append(e)
        if len(I) == RANK: break
    assert len(I) == RANK
    OUT = [e for e in S.SAT if e not in I]; NOUT = len(OUT)
    AI = [ROWS[e] for e in I]
    aug = [AI[i] + [1 if k == i else 0 for k in range(RANK)] for i in range(RANK)]
    R, piv = gauss(aug, NK)
    C = [[0]*NK for _ in range(RANK)]
    for i, c in enumerate(piv):
        for j in range(RANK): C[j][c] = R[i][NK + j]
    for i, e in enumerate(I):
        for j in range(RANK):
            v = sum(a*b for a, b in zip(ROWS[e], C[j])) % P
            assert v == (1 if i == j else 0)
    G = [[sum(a*b for a, b in zip(ROWS[e], C[j])) % P for e in OUT] for j in range(RANK)]
    before = len(found); degen = 0
    def record(Tidx, y):
        nz = [k for k in range(NOUT) if sum(y[a]*G[j][k] for a, j in enumerate(Tidx)) % P]
        if len(nz) <= 6 - len(Tidx):
            found.add(frozenset([I[j] for j in Tidx] + [OUT[k] for k in nz])); return 1
        return 0
    for j in range(RANK): record((j,), [1])
    for i, j in itertools.combinations(range(RANK), 2):
        ratios = Counter(); always = 0
        for k in range(NOUT):
            a, b = G[i][k], G[j][k]
            if a == 0 and b == 0: always += 1
            elif b: ratios[(-a)*pow(b, P-2, P) % P] += 1
        for lam, m in ratios.most_common(4):
            if always + m >= NOUT - 4: record((i, j), [1, lam])
    print('  [set %d] s<=2 done: %d supports, %.0fs' % (trial, len(found), time.time()-t0), flush=True)
    for s in range(3, 7):
        WIN = 8 if s <= 4 else 6
        for Tidx in itertools.combinations(range(RANK), s):
            win = []
            for k in range(NOUT):
                if any(G[j][k] for j in Tidx):
                    win.append(k)
                    if len(win) == WIN: break
            for Z in itertools.combinations(win, s-1):
                M = [[G[j][k] for j in Tidx] for k in Z]
                Rr, pv = gauss(M, s)
                free = [c for c in range(s) if c not in pv]
                if len(free) != 1: degen += 1; continue
                f = free[0]; y = [0]*s; y[f] = 1
                for i2, c in enumerate(pv): y[c] = (-Rr[i2][f]) % P
                record(Tidx, y)
        print('  [set %d] s=%d done: %d supports, %.0fs' % (trial, s, len(found), time.time()-t0), flush=True)
    degen_tot += degen
    yields.append(len(found) - before)
    print('[information set %d] NEW supports: %d  (total %d, degenerate skipped %d)'
          % (trial, len(found) - before, len(found), degen), flush=True)
print('\nmarginal yield per information set:', yields)
print('degenerate (skipped) subsets total: %d' % degen_tot)
json.dump({'supports': [sorted(map(str, T)) for T in sorted(found, key=lambda t: (len(t), sorted(map(str, t))))],
           'yields': yields, 'degenerate': degen_tot},
          open('w_cocirc2_raw.json', 'w'), indent=1)
print('raw candidate supports: %d  elapsed %.0fs' % (len(found), time.time()-t0))
