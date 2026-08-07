#!/usr/bin/env python3
"""agent AF, step 15: full anatomy of one block; and the global freedom census."""
import sys, os, pickle
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, pp, expand, varsof, shape_of, Pval
from af1_parse import is_const
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
B = pickle.load(open(os.path.join(HERE, 'af_blocks.pkl'), 'rb'))
LAM = pickle.load(open(os.path.join(HERE, 'af_lam.pkl'), 'rb'))
conds = C['conds']; info = M['info']; blk = M['blk']; live = M['live']; pure = LAM['pure']

def raw(n, d=0):
    if n[0] == 'c':
        return str(n[1]) if abs(n[1]) < 10**14 else 'P' if n[1] == Pval else 'BIG'
    if n[0] == 'v':
        return 'x%d' % find(n[1])
    if n[0] == 'neg':
        return '-' + raw(n[1])
    return '(%s %s %s)' % (raw(n[1]), n[0], raw(n[2]))

# --- collect the diff conditions and index them by (X,Y)
diffidx = {}
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    if info[i]['cls'] == 'diff':
        diffidx[(info[i]['X'], info[i]['Y'])] = (i, c)

# --- for each block, the 3 congruences: residual B = a1*(X1-Y1) + a2*(X2-Y2)
blkdata = {}
for i, (R, c, Ex, aid, uc) in enumerate(conds):
    d = info[i]
    if d['cls'] != 'cong':
        continue
    dl = defc[d['other']]
    rhs = dl[0][1]
    assert rhs[0] == '+'
    def term(t):
        if t[0] == '*':
            k = is_const(t[1])
            if k is not None:
                return (k, t[2])
            k = is_const(t[2])
            if k is not None:
                return (k, t[1])
        if t[0] == 'neg':
            k, x = term(t[1]); return (-k, x)
        return (1, t)
    (a1, D1), (a2, D2) = term(rhs[1]), term(rhs[2])
    def pair(D):
        if D[0] == 'v':
            e = expand(D, 1)
            if e[0] == '-' and e[1][0] == 'v' and e[2][0] == 'v':
                return (find(e[1][1]), find(e[2][1]))
            return ('w', find(D[1]))
        if D[0] == '-' and D[1][0] == 'v' and D[2][0] == 'v':
            return (find(D[1][1]), find(D[2][1]))
        return ('?', shape_of(D))
    blkdata.setdefault(d['gate'], []).append((c, a1, pair(D1), a2, pair(D2)))
print('blocks with 3 congruence rows: %d' % sum(1 for v in blkdata.values() if len(v) == 3))

# do all 3 rows of a block share the SAME two differences (N1,N2)?
same = 0; nd = Counter()
for g, rows in blkdata.items():
    ds = set()
    for (c, a1, D1, a2, D2) in rows:
        ds.add(D1); ds.add(D2)
    nd[len(ds)] += 1
    if len(ds) == 2:
        same += 1
print('blocks whose 3 rows use exactly 2 distinct differences: %d ; histogram %s' % (same, dict(nd)))

# are those differences exactly the ungated `diff` conditions?
hit = 0; tot = 0
for g, rows in blkdata.items():
    ds = set()
    for (c, a1, D1, a2, D2) in rows:
        ds.add(D1); ds.add(D2)
    for D in ds:
        tot += 1
        if D in diffidx:
            hit += 1
print('block differences that are ALSO an ungated  c*P | (X-Y)  condition: %d / %d' % (hit, tot))

# rank of the 3x2 coefficient matrix, over Q and the size of minors
import itertools
bad = 0; mx = 0
for g, rows in blkdata.items():
    ds = []
    for (c, a1, D1, a2, D2) in rows:
        for D in (D1, D2):
            if D not in ds:
                ds.append(D)
    if len(ds) != 2:
        continue
    Mx = []
    for (c, a1, D1, a2, D2) in rows:
        r = [0, 0]
        r[ds.index(D1)] += a1; r[ds.index(D2)] += a2
        Mx.append(r)
    mins = [Mx[i][0]*Mx[j][1] - Mx[j][0]*Mx[i][1] for i, j in itertools.combinations(range(3), 2)]
    if any(m == 0 for m in mins):
        bad += 1
    mx = max(mx, max(abs(m) for m in mins))
print('blocks with a vanishing 2x2 minor: %d ;  max |minor| = %d  (< P: %s)' % (bad, mx, mx < Pval))

# ---------- one live block, in full ----------
g = sorted(live, key=lambda x: len(pure[x][0]) + len(pure[x][1]))[0]
I, J = pure[g]
print('\n===== smallest merge block: gate x%d,  I=%s  J=%s' % (g, sorted(I), sorted(J)))
for (c, a1, D1, a2, D2) in blkdata[g]:
    print('   c = %-12d   %d*(x%d - x%d)  +  %d*(x%d - x%d)' % (c, a1, D1[0], D1[1], a2, D2[0], D2[1]))
for i, (R, cc, Ex, aid, uc) in enumerate(conds):
    if info[i]['cls'] == 'offpin' and info[i]['gate'] == g:
        print('   off-pin  c = %-12d  on free wire x%d' % (cc, info[i]['other']))
ds = set()
for (c, a1, D1, a2, D2) in blkdata[g]:
    ds.add(D1); ds.add(D2)
for D in ds:
    for w in D:
        print('   x%d = %s' % (w, raw(expand(('v', w), 4))[:260]))
    if D in diffidx:
        print('     ungated condition  %d*P | (x%d - x%d)' % (diffidx[D][1], D[0], D[1]))
pickle.dump({'blkdata': blkdata, 'diffidx': diffidx}, open(os.path.join(HERE, 'af_blkdata.pkl'), 'wb'), 2)
