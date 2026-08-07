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
    if Ex[0] == '*' and Ex[1][0] == 'v' and Ex[2][0] == '-' \
       and Ex[2][1][0] == 'v' and is_const(Ex[2][2]) is not None:
        leafpin.append((find(Ex[1][1]), find(Ex[2][1][1]), is_const(Ex[2][2]), c, R))
sel = sorted(set(t[0] for t in leafpin))
print('leaf-pin conditions: %d   distinct selectors: %d' % (len(leafpin), len(sel)))

# ---------- 2. booleanity: any atom equivalent to v^2 - v ----------
boolvar = set()
def bool_from(a):
    if a[0] != '-':
        return None
    X, Y = a[1], a[2]
    def sq(n):
        if n[0] == '*' and n[1][0] == 'v' and n[2][0] == 'v' and find(n[1][1]) == find(n[2][1]):
            return find(n[1][1])
        return None
    def lin(n):
        if n[0] == 'v':
            return (1, find(n[1]))
        if n[0] == '*':
            k = is_const(n[1])
            if k is not None and n[2][0] == 'v':
                return (k, find(n[2][1]))
            k = is_const(n[2])
            if k is not None and n[1][0] == 'v':
                return (k, find(n[1][1]))
        return None
    s, l = sq(X), lin(Y)
    if s is not None and l is not None and l[1] == s and abs(l[0]) == 1:
        return s
    s, l = sq(Y), lin(X)
    if s is not None and l is not None and l[1] == s and abs(l[0]) == 1:
        return s
    return None
for a in atoms:
    b = bool_from(a)
    if b is not None:
        boolvar.add(b)
    # v*(v-1)
    if a[0] == '*' and a[1][0] == 'v' and a[2][0] == '-' and a[2][1][0] == 'v' \
       and find(a[2][1][1]) == find(a[1][1]) and is_const(a[2][2]) == 1:
        boolvar.add(find(a[1][1]))
    if a[0] == '*' and a[2][0] == 'v' and a[1][0] == '-' and a[1][1][0] == 'v' \
       and find(a[1][1][1]) == find(a[2][1]) and is_const(a[1][2]) == 1:
        boolvar.add(find(a[2][1]))
print('variables with an explicit booleanity atom: %d' % len(boolvar))
print('  selectors that are boolean: %d / %d' % (len(set(sel) & boolvar), len(sel)))

# ---------- 3. gate definitions ----------
kind = {}
for r, dl in defc.items():
    for aid, rhs in dl:
        k = None
        if rhs[0] == '*' and rhs[1][0] == 'v' and rhs[2][0] == 'v':
            k = ('AND', find(rhs[1][1]), find(rhs[2][1]))
        elif rhs[0] == '-' and rhs[1][0] == '+' and rhs[2][0] == '*' \
             and rhs[1][1][0] == 'v' and rhs[1][2][0] == 'v' \
             and rhs[2][1][0] == 'v' and rhs[2][2][0] == 'v' \
             and {find(rhs[1][1][1]), find(rhs[1][2][1])} == {find(rhs[2][1][1]), find(rhs[2][2][1])}:
            k = ('OR', find(rhs[1][1][1]), find(rhs[1][2][1]))
        elif rhs[0] == '-' and is_const(rhs[1]) == 1 and rhs[2][0] == 'v':
            k = ('NOT', find(rhs[2][1]))
        elif rhs[0] == '+' and rhs[1][0] == 'v' and rhs[2][0] == 'v':
            k = ('SUM', find(rhs[1][1]), find(rhs[2][1]))
        elif rhs[0] == 'c':
            k = ('CONST', rhs[1])
        if k is not None:
            kind.setdefault(r, k)

selset = set(sel)
known = dict((s, ('LEAF', s)) for s in selset)
for r, x in kind.items():
    if x[0] == 'CONST' and x[1] in (0, 1):
        known[r] = x
changed = True; rounds = 0
while changed and rounds < 60:
    changed = False; rounds += 1
    for r, x in kind.items():
        if r in known:
            continue
        if x[0] == 'NOT':
            if x[1] in known:
                known[r] = x; changed = True
        elif x[0] in ('AND', 'OR', 'SUM'):
            if x[1] in known and x[2] in known:
                known[r] = x; changed = True
print('boolean-cone nodes reachable from selectors: %d (rounds %d)' % (len(known), rounds))
print('  kinds:', dict(Counter(v[0] for v in known.values())))

# ---------- 4. support of each boolean node (set of selectors it depends on) ----------
import functools
order = []
seen = set()
def topo(r):
    st = [(r, 0)]
    while st:
        n, ph = st.pop()
        if ph == 0:
            if n in seen:
                continue
            seen.add(n)
            st.append((n, 1))
            x = known.get(n)
            if x and x[0] != 'LEAF' and x[0] != 'CONST':
                for c in x[1:]:
                    if c in known:
                        st.append((c, 0))
        else:
            order.append(n)
for r in known:
    topo(r)
supp = {}
for r in order:
    x = known[r]
    if x[0] == 'LEAF':
        supp[r] = frozenset([x[1]])
    elif x[0] == 'CONST':
        supp[r] = frozenset()
    elif x[0] == 'NOT':
        supp[r] = supp[x[1]]
    else:
        supp[r] = supp[x[1]] | supp[x[2]]
sizes = Counter(len(s) for s in supp.values())
print('  distinct supports among boolean nodes: %d' % len(set(supp.values())))
print('  support size histogram (top):', sorted(sizes.items())[:12], '...max', max(sizes))

pickle.dump({'sel': sel, 'leafpin': leafpin, 'boolvar': boolvar,
             'kind': kind, 'known': known, 'supp': supp},
            open(os.path.join(HERE, 'af_tree.pkl'), 'wb'), 2)
