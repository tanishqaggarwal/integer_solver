#!/usr/bin/env python3
"""agent AF, step 18: the per-block freedom ledger, exhaustively over all 383 blocks.

Checks, for every block:
  (a) its two off-pin wires are FREE and are used only inside its own block cluster;
  (b) the 3x2 congruence coefficient matrix has rank 2 over Q, with all minors < P;
  (c) at most one congruence row and at most one off-pin row carries c > 1;
  (d) gcd(c, alpha) and gcd(c, beta) for the c>1 rows;
  (e) what the 766 ungated `c*P | X - Y` conditions are.
"""
import sys, os, pickle, itertools
from math import gcd
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, pp, expand, varsof, shape_of, Pval
from af1_parse import is_const
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
R = pickle.load(open(os.path.join(HERE, 'af_rows.pkl'), 'rb'))
conds = C['conds']; info = M['info']; rows = R['rows']

v2a = defaultdict(list)
for aid, a in enumerate(atoms):
    for v in varsof(a, set()):
        v2a[find(v)].append(aid)

# ---------- (a) off-pin wires ----------
opw = defaultdict(list)
for i, (Rw, c, Ex, aid, uc) in enumerate(conds):
    if info[i]['cls'] == 'offpin':
        opw[info[i]['gate']].append((info[i]['other'], c))
free = sum(1 for g in opw for (w, c) in opw[g] if not defc.get(w))
print('(a) off-pin wires: %d ; with zero definitions (FREE): %d' %
      (sum(len(v) for v in opw.values()), free))
# is the wire ever multiplied by its own gate L (the mux)?
muxed = 0
for g, lst in opw.items():
    for (w, c) in lst:
        hit = False
        for aid in v2a[w]:
            a = atoms[aid]
            if a[0] == '-' and a[2][0] == '*' and a[2][1][0] == 'v' and a[2][2][0] == 'v' \
               and {find(a[2][1][1]), find(a[2][2][1])} == {g, w}:
                hit = True
        muxed += hit
print('    off-pin wires that appear in a mux atom  V - (L * wire): %d' % muxed)
short = Counter()
for g, lst in opw.items():
    for (w, c) in lst:
        n = sum(1 for aid in v2a[w] if len(str(atoms[aid])) < 400)
        short[n] += 1
print('    #SHORT atoms touching each off-pin wire:', dict(sorted(short.items())))

# ---------- (b),(c),(d) congruence matrix ----------
# canonical id of a difference wire: expand fully and hash the leaf variables
def diffid(D):
    if D[0] == 'v':
        e = expand(D, 1)
        if e[0] == '-' and e[1][0] == 'v' and e[2][0] == 'v':
            return ('d', find(e[1][1]), find(e[2][1]))
        return ('w', find(D[1]))
    if D[0] == '-' and D[1][0] == 'v' and D[2][0] == 'v':
        return ('d', find(D[1][1]), find(D[2][1]))
    return ('?', shape_of(D))

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

# Each row's two differences must be identified ACROSS rows of the same block.
# Do it by the pair of *ultimate* wires: expand the difference operands one more level.
def deep(v, d=3):
    return expand(('v', v), d)

blkrows = defaultdict(list)
for i, (Rw, c, Ex, aid, uc) in enumerate(conds):
    if info[i]['cls'] != 'cong':
        continue
    rhs = defc[info[i]['other']][0][1]
    (a1, D1), (a2, D2) = term(rhs[1]), term(rhs[2])
    blkrows[info[i]['gate']].append((c, a1, diffid(D1), a2, diffid(D2)))

# identify differences across rows by their fully-expanded printed form
def canon(did):
    if did[0] == 'd':
        return (shape_of(deep(did[1], 6)), shape_of(deep(did[2], 6)),
                tuple(sorted(varsof(deep(did[1], 6), set()) | varsof(deep(did[2], 6), set()))))
    return did

rk = Counter(); mx = 0; badminor = 0
for g, rs in blkrows.items():
    ids = []
    for (c, a1, d1, a2, d2) in rs:
        for d in (d1, d2):
            k = canon(d)
            if k not in ids:
                ids.append(k)
    if len(ids) != 2:
        rk[len(ids)] += 1
        continue
    rk[2] += 1
    Mx = []
    for (c, a1, d1, a2, d2) in rs:
        r = [0, 0]
        r[ids.index(canon(d1))] += a1
        r[ids.index(canon(d2))] += a2
        Mx.append(r)
    mins = [Mx[i][0]*Mx[j][1] - Mx[j][0]*Mx[i][1] for i, j in itertools.combinations(range(3), 2)]
    if any(m == 0 for m in mins):
        badminor += 1
    mx = max(mx, max(abs(m) for m in mins))
print('\n(b) #distinct differences per block:', dict(rk))
print('    blocks with a vanishing 2x2 minor: %d ; max |minor| = %d ; < P: %s' % (badminor, mx, mx < Pval))

print('\n(c) per block: #cong rows c>1 :', dict(sorted(Counter(sum(1 for r in rs if r[0] > 1) for rs in blkrows.values()).items())))
print('    per block: #offpin rows c>1:', dict(sorted(Counter(sum(1 for (w, c) in lst if c > 1) for lst in opw.values()).items())))

print('\n(d) for the 288 c>1 congruence rows:')
ga = Counter(); gb = Counter()
for (g, c, a1, a2, D1, D2) in rows:
    if c > 1:
        ga[gcd(c, a1)] += 1; gb[gcd(c, a2)] += 1
print('    gcd(c,alpha):', ga.most_common(5))
print('    gcd(c,beta) :', gb.most_common(5))
print('    rows where gcd(c,alpha)=1 OR gcd(c,beta)=1: %d / 288' %
      sum(1 for (g, c, a1, a2, D1, D2) in rows if c > 1 and (gcd(c, a1) == 1 or gcd(c, a2) == 1)))

# ---------- (e) the 766 ungated differences ----------
print('\n(e) the 766 ungated  c*P | (X - Y)  conditions:')
def raw(n):
    if n[0] == 'c':
        return str(n[1]) if abs(n[1]) < 10**14 else ('P' if n[1] == Pval else 'BIG')
    if n[0] == 'v':
        return 'x%d' % find(n[1])
    if n[0] == 'neg':
        return '-' + raw(n[1])
    return '(%s %s %s)' % (raw(n[1]), n[0], raw(n[2]))
dfree = Counter()
n = 0
for i, (Rw, c, Ex, aid, uc) in enumerate(conds):
    if info[i]['cls'] != 'diff':
        continue
    X, Y = info[i]['X'], info[i]['Y']
    dfree[(len(defc.get(X, [])), len(defc.get(Y, [])))] += 1
    if n < 5:
        print('    c=%-10d  x%d - x%d   with  x%d = %s ,  x%d = %s'
              % (c, X, Y, X, raw(expand(('v', X), 2))[:90], Y, raw(expand(('v', Y), 2))[:90]))
        n += 1
print('    (#defs of X, #defs of Y) census:', dict(dfree))
