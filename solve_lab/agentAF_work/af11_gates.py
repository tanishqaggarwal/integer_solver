#!/usr/bin/env python3
"""agent AF, step 11: expand the 383 liveness gates as boolean functions of the 256 selectors."""
import sys, os, pickle
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af6_expand import (atoms, defs, defc, val, lift, Pval, find, pp, expand, varsof, shape_of)
from af1_parse import is_const
B = pickle.load(open(os.path.join(HERE, 'af_blocks.pkl'), 'rb'))
T = pickle.load(open(os.path.join(HERE, 'af_tree.pkl'), 'rb'))
gates = B['gates']; sel = T['sel']; selset = set(sel)
selidx = {s: i for i, s in enumerate(sel)}

# node types: ('leaf',i) ('and',a,b) ('or',a,b) ('not',a) ('const',k) ('sum',a,b) ('opq',r)
memo = {}
STACK = set()
OPQ = ('opq',)

def mk_and(a, b):
    if a[0] == 'const':
        return b if a[1] == 1 else ('const', 0)
    if b[0] == 'const':
        return a if b[1] == 1 else ('const', 0)
    if a == b:
        return a
    return ('and',) + tuple(sorted([a, b], key=repr))

def mk_or(a, b):
    if a[0] == 'const':
        return ('const', 1) if a[1] == 1 else b
    if b[0] == 'const':
        return ('const', 1) if b[1] == 1 else a
    if a == b:
        return a
    return ('or',) + tuple(sorted([a, b], key=repr))

def bexp(r, depth=0):
    r = find(r)
    if r in memo:
        return memo[r]
    if r in selset:
        memo[r] = ('leaf', selidx[r]); return memo[r]
    if r in STACK or depth > 600:
        return ('opq', r)
    dl = defc.get(r)
    if not dl:
        return ('opq', r)
    STACK.add(r)
    res = ('opq', r)
    for aid, rhs in dl:
        v = bnode(rhs, depth+1)
        if v is not None and v[0] != 'opq':
            res = v; break
    STACK.discard(r)
    memo[r] = res
    return res

def bnode(n, depth):
    t = n[0]
    if t == 'c':
        return ('const', n[1])
    if t == 'v':
        return bexp(n[1], depth)
    if t == 'neg':
        a = bnode(n[1], depth)
        if a is not None and a[0] == 'const':
            return ('const', -a[1])
        return None
    a = bnode(n[1], depth); b = bnode(n[2], depth)
    if a is None or b is None:
        return None
    if t == '*':
        if a[0] == 'opq' or b[0] == 'opq':
            if a[0] == 'const' and a[1] == 0:
                return ('const', 0)
            if b[0] == 'const' and b[1] == 0:
                return ('const', 0)
            return None
        return mk_and(a, b)
    if t == '+':
        if a[0] == 'const' and b[0] == 'const':
            return ('const', a[1] + b[1])
        if a[0] == 'const' and a[1] == 0:
            return b
        if b[0] == 'const' and b[1] == 0:
            return a
        if a[0] == 'opq' or b[0] == 'opq':
            return None
        return ('sum',) + tuple(sorted([a, b], key=repr))
    if t == '-':
        if a[0] == 'const' and b[0] == 'const':
            return ('const', a[1] - b[1])
        if a[0] == 'const' and a[1] == 1 and b[0] != 'opq':
            return b[1] if b[0] == 'not' else ('not', b)
        # sum(x,y) - and(x,y)  ->  or(x,y)
        if a[0] == 'sum' and b[0] == 'and' and set(a[1:]) == set(b[1:]):
            return mk_or(a[1], a[2])
        if b[0] == 'const' and b[1] == 0:
            return a
        return None
    return None

out = {}
for L in gates:
    out[L] = bexp(L)
print('gate expansion kinds:', dict(Counter(v[0] for v in out.values())))

def support(n, acc):
    if n[0] == 'leaf':
        acc.add(n[1])
    elif n[0] in ('and', 'or', 'sum'):
        support(n[1], acc); support(n[2], acc)
    elif n[0] == 'not':
        support(n[1], acc)
    return acc

def depthof(n, memo={}):
    if n[0] in ('leaf', 'const', 'opq'):
        return 0
    if n[0] == 'not':
        return 1 + depthof(n[1])
    return 1 + max(depthof(n[1]), depthof(n[2]))

good = {L: f for L, f in out.items() if f[0] != 'opq'}
print('fully expanded gates: %d / %d' % (len(good), len(out)))
print('const gates:', dict(Counter(f[1] for f in out.values() if f[0] == 'const')))
print('support-size histogram:', dict(sorted(Counter(len(support(f, set())) for f in good.values()).items())))
print('depth histogram:', dict(sorted(Counter(depthof(f) for f in good.values()).items())))
tc = Counter()
def census(n, seen):
    if id(n) in seen:
        return
    seen.add(id(n)); tc[n[0]] += 1
    if n[0] in ('and', 'or', 'sum'):
        census(n[1], seen); census(n[2], seen)
    elif n[0] == 'not':
        census(n[1], seen)
for f in good.values():
    census(f, set())
print('operator census inside gate formulas:', dict(tc))

pickle.dump({'gatefn': out, 'selidx': selidx, 'sel': sel},
            open(os.path.join(HERE, 'af_gates.pkl'), 'wb'), 2)
