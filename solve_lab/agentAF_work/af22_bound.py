#!/usr/bin/env python3
"""agent AF, step 22: which variables are BOUNDED?  (a growth argument needs bounded variables)"""
import sys, os, pickle
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import atoms, defc, val, find, pp, expand, varsof, shape_of, Pval
from af1_parse import is_const
C = pickle.load(open(os.path.join(HERE, 'af_cond.pkl'), 'rb'))
M = pickle.load(open(os.path.join(HERE, 'af_map.pkl'), 'rb'))
T = pickle.load(open(os.path.join(HERE, 'af_tree.pkl'), 'rb'))
LAM = pickle.load(open(os.path.join(HERE, 'af_lam.pkl'), 'rb'))
conds = C['conds']; info = M['info']; sel = T['sel']; pure = LAM['pure']

# --- booleanity, exhaustive over every atom shape that forces v(v-1)=0 ---
boolv = set()
def lin(n):
    if n[0] == 'v':
        return (1, find(n[1]))
    if n[0] == 'neg':
        r = lin(n[1]); return None if r is None else (-r[0], r[1])
    if n[0] == '*':
        k = is_const(n[1])
        if k is not None and n[2][0] == 'v':
            return (k, find(n[2][1]))
        k = is_const(n[2])
        if k is not None and n[1][0] == 'v':
            return (k, find(n[1][1]))
    return None
def sq(n):
    if n[0] == '*' and n[1][0] == 'v' and n[2][0] == 'v' and find(n[1][1]) == find(n[2][1]):
        return find(n[1][1])
    if n[0] == '*':
        k = is_const(n[1])
        if k is not None:
            s = sq(n[2]); return s
        k = is_const(n[2])
        if k is not None:
            return sq(n[1])
    return None
for a in atoms:
    if a[0] == '-':
        for X, Y in ((a[1], a[2]), (a[2], a[1])):
            s, l = sq(X), lin(Y)
            if s is not None and l is not None and l[1] == s:
                boolv.add(s)
    if a[0] == '*':
        for X, Y in ((a[1], a[2]), (a[2], a[1])):
            if X[0] == 'v' and Y[0] == '-':
                if Y[1][0] == 'v' and find(Y[1][1]) == find(X[1]) and is_const(Y[2]) == 1:
                    boolv.add(find(X[1]))
                if is_const(Y[1]) == 1 and Y[2][0] == 'v' and find(Y[2][1]) == find(X[1]):
                    boolv.add(find(X[1]))
print('variables forced boolean by an explicit atom: %d' % len(boolv))
print('  of the 256 selectors: %d' % len(set(find(s) for s in sel) & boolv))
opw = set(info[i]['other'] for i in range(len(conds)) if info[i]['cls'] == 'offpin')
print('  of the 766 free chord-output wires: %d' % len(opw & boolv))
lw = set()
for i in range(len(conds)):
    if info[i]['cls'] == 'pin':
        E = conds[i][2]
        lw.add(find(E[2][1][1]))
print('  of the %d leaf coordinate wires: %d' % (len(lw), len(lw & boolv)))

# --- support distribution of the 310 merge-block c>1 conditions ---
merge = set(pure)
szs = []
for i in range(len(conds)):
    if conds[i][1] == 1:
        continue
    d = info[i]
    if d['cls'] in ('cong', 'offpin') and d['gate'] in merge:
        I, J = pure[d['gate']]
        szs.append(len(I) + len(J))
szs.sort()
print('\nsupport sizes of the 310 gated c>1 conditions on merge blocks:')
print('  ', Counter(szs).most_common(10))
for th in (1, 2, 4, 8, 16, 32, 64, 128, 255):
    print('   support > %3d : %3d conditions' % (th, sum(1 for s in szs if s > th)))

# --- how many conditions can even see |S| at all? ---
print('\nconditions whose selector support is the WHOLE leaf set (256): %d'
      % sum(1 for s in szs if s == 256))
