#!/usr/bin/env python3
"""agent AF, step 9: selectors, the boolean liveness circuit, the 383 blocks, the tree."""
import sys, os, pickle, json
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import (atoms, defs, defc, val, lift, Pval, find, pp, expand, varsof, shape_of)
from af1_parse import is_const
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
conds = C['conds']

# ---------- 1. leaf pins -> selectors ----------
leafpin = []
for (R, c, Ex, aid, uc) in conds:
    if Ex[0] == '*' and Ex[1][0] == 'v' and Ex[2][0] == '-' and Ex[2][0] == '-' \
       and Ex[2][1][0] == 'v' and is_const(Ex[2][2]) is not None:
        leafpin.append((find(Ex[1][1]), find(Ex[2][1][1]), is_const(Ex[2][2]), c, R))
sel = sorted(set(t[0] for t in leafpin))
print('leaf-pin conditions: %d   distinct selectors: %d' % (len(leafpin), len(sel)))
byc = Counter(t[3] for t in leafpin)
print('  with c>1: %d ; with c==1: %d' % (sum(k for c, k in byc.items() if c > 1), byc.get(1, 0)))
persel = Counter(t[0] for t in leafpin)
print('  pins per selector:', dict(Counter(persel.values())))

# ---------- 2. booleanity ----------
boolvar = set()
for a in atoms:
    # (V*V) - V   /  V - (V*V)  / V*(V-1)
    if a[0] == '-' and a[1][0] == '*' and a[1][1] == a[1][2] and a[2][0] == 'v' \
       and find(a[1][1][1]) == find(a[2][1]):
        boolvar.add(find(a[2][1]))
    if a[0] == '-' and a[2][0] == '*' and a[2][1] == a[2][2] and a[1][0] == 'v' \
       and find(a[2][1][1]) == find(a[1][1]):
        boolvar.add(find(a[1][1]))
    if a[0] == '*' and a[1][0] == 'v' and a[2][0] == '-' and a[2][1][0] == 'v' \
       and find(a[2][1][1]) == find(a[1][1]) and is_const(a[2][2]) == 1:
        boolvar.add(find(a[1][1]))
print('variables with an explicit booleanity atom: %d' % len(boolvar))
print('  selectors that are boolean: %d / %d' % (len(set(sel) & boolvar), len(sel)))

# ---------- 3. boolean circuit above the selectors ----------
# defs of the form  v = a*b  (AND),  v = a+b (SUM),  v = 1 - a (NOT),  v = K
kind = {}
for r, dl in defc.items():
    if len(dl) != 1:
        continue
    rhs = dl[0][1]
    if rhs[0] == '*' and rhs[1][0] == 'v' and rhs[2][0] == 'v':
        kind[r] = ('AND', find(rhs[1][1]), find(rhs[2][1]))
    elif rhs[0] == '+' and rhs[1][0] == 'v' and rhs[2][0] == 'v':
        kind[r] = ('SUM', find(rhs[1][1]), find(rhs[2][1]))
    elif rhs[0] == '-' and is_const(rhs[1]) == 1 and rhs[2][0] == 'v':
        kind[r] = ('NOT', find(rhs[2][1]))
    elif rhs[0] == '-' and rhs[1][0] == 'v' and rhs[2][0] == 'v':
        kind[r] = ('SUB', find(rhs[1][1]), find(rhs[2][1]))

# forward closure: which nodes are computable from selectors alone through AND/SUM/NOT/SUB
selset = set(sel)
known = dict((s, ('LEAF', s)) for s in selset)
for r, x in kind.items():
    pass
changed = True
rounds = 0
while changed and rounds < 40:
    changed = False; rounds += 1
    for r, x in kind.items():
        if r in known:
            continue
        if x[0] == 'NOT':
            if x[1] in known:
                known[r] = x; changed = True
        else:
            if x[1] in known and x[2] in known:
                known[r] = x; changed = True
print('boolean-cone nodes reachable from selectors: %d (rounds %d)' % (len(known), rounds))
print('  kinds:', dict(Counter(v[0] for v in known.values())))

pickle.dump({'sel': sel, 'leafpin': leafpin, 'boolvar': boolvar,
             'kind': kind, 'known': known},
            open(os.path.join(HERE, 'af_tree.pkl'), 'wb'), 2)
