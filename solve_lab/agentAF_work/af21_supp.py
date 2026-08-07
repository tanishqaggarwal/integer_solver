#!/usr/bin/env python3
"""agent AF, step 21: exact selector-support of all 3707 conditions; the 927 table."""
import sys, os, pickle
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, pp, expand, varsof, shape_of, Pval
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
LAM = pickle.load(open(os.path.join(HERE, 'af_lam.pkl'), 'rb'))
T = pickle.load(open(os.path.join(HERE, 'af_tree.pkl'), 'rb'))
conds = C['conds']; info = M['info']; pure = LAM['pure']; sel = T['sel']
selidx = {s: i for i, s in enumerate(sel)}

# ---- (1) the 180 alias-identical diffs, re-checked ----
copyset = set()
for a in atoms:
    if a[0] == '-' and a[1][0] == 'v' and a[2][0] == 'v':
        copyset.add(frozenset((a[1][1], a[2][1])))
tot = ok = 0
for i in range(len(conds)):
    if info[i]['cls'] != 'diff' or info[i]['X'] != info[i]['Y']:
        continue
    tot += 1
    E = conds[i][2]
    if frozenset((E[1][1], E[2][1])) in copyset:
        ok += 1
print('(1) diff conditions with X and Y in the same alias class: %d ; each backed by an explicit copy atom: %d' % (tot, ok))
nfreeX = sum(1 for i in range(len(conds))
             if info[i]['cls'] == 'diff' and info[i]['X'] != info[i]['Y']
             and not defc.get(info[i]['X']))
print('    remaining diff conditions whose LHS wire X is FREE (0 definitions): %d' % nfreeX)
print('    -> vacuous %d + free-LHS %d = %d of 766' % (tot, nfreeX, tot + nfreeX))

# ---- (2) exact selector support of every condition ----
# forward selector-support closure over the whole definition DAG
SUP = {}
def supp(v, depth=0):
    v = find(v)
    if v in SUP:
        return SUP[v]
    if v in selidx:
        SUP[v] = frozenset([selidx[v]]); return SUP[v]
    if v in val:
        SUP[v] = frozenset(); return SUP[v]
    dl = defc.get(v)
    if not dl or depth > 400:
        SUP[v] = frozenset(); return SUP[v]
    SUP[v] = frozenset()
    acc = set()
    for aid, rhs in dl:
        for w in varsof(rhs, set()):
            acc |= supp(w, depth + 1)
    SUP[v] = frozenset(acc)
    return SUP[v]

sys.setrecursionlimit(100000)
for v in list(defc):
    supp(v)

merge = set(pure)
tbl = defaultdict(list)
for i, (Rw, c, Ex, aid, uc) in enumerate(conds):
    if c == 1:
        continue
    d = info[i]; cls = d['cls']
    if cls == 'pin':
        s = frozenset([selidx[d['selector']]]) if d['selector'] in selidx else frozenset()
        tag = 'leaf pin'
    elif cls in ('cong', 'offpin'):
        g = d['gate']
        if g in merge:
            I, J = pure[g]
            s = I | J
            tag = '%s @ merge block' % cls
        else:
            s = frozenset()
            tag = '%s @ constant-0-gate block' % cls
    elif cls == 'diff':
        if d['X'] == d['Y']:
            s = frozenset(); tag = 'diff (identically 0)'
        else:
            s = supp(d['X']) | supp(d['Y']); tag = 'diff (free LHS)'
    else:
        s = frozenset(); tag = cls
    tbl[tag].append(len(s))
print('\n(2) selector support of the 927 conditions with c > 1:')
tot927 = 0
for k in sorted(tbl):
    v = tbl[k]; tot927 += len(v)
    h = Counter(v)
    print('    %-34s n=%4d   support sizes: min %3d  median %3d  max %3d'
          % (k, len(v), min(v), sorted(v)[len(v)//2], max(v)))
print('    total %d' % tot927)
big = [(k, sorted(v)[-3:]) for k, v in tbl.items() if max(v) > 8]
print('    families with support > 8:', big)
