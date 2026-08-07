#!/usr/bin/env python3
"""agent AF, step 6: symbolic expansion of the lift atoms down to free variables."""
import sys, os, pickle
from collections import Counter, defaultdict
sys.setrecursionlimit(200000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af2_atoms import shape_of, varsof

D = pickle.load(open(os.path.join(HERE, 'af_atoms.pkl'), 'rb'))
E = pickle.load(open(os.path.join(HERE, 'af_defs.pkl'), 'rb'))
F = pickle.load(open(os.path.join(HERE, 'af_P.pkl'), 'rb'))
G = pickle.load(open(os.path.join(HERE, 'af_lift.pkl'), 'rb'))
atoms = D['atoms']; defs = E['defs']; Pval = F['Pval']
parent = F['parent']; cval = F['cval']; val = G['val']; lift = G['lift']

def find(x):
    parent.setdefault(x, x)
    r = x
    while parent[r] != r:
        r = parent[r]
    while parent[x] != r:
        parent[x], x = r, parent[x]
    return r

# canonical definition per class: prefer non-copy definitions
defc = {}
for v, lst in defs.items():
    r = find(v)
    for aid, rhs in lst:
        if rhs[0] == 'v' and find(rhs[1]) == r:
            continue
        defc.setdefault(r, []).append((aid, rhs))

def ndefs(v):
    return len(defc.get(find(v), ()))

def pp(n, depth=0):
    if n[0] == 'c':
        return str(n[1])
    if n[0] == 'v':
        r = find(n[1])
        if r in val:
            return 'K[%d]' % val[r]
        return 'x%d' % r
    if n[0] == 'neg':
        return '-' + pp(n[1], depth+1)
    return '(%s %s %s)' % (pp(n[1], depth+1), n[0], pp(n[2], depth+1))

def expand(n, depth, stop=frozenset()):
    """unfold variable definitions to given depth"""
    if depth <= 0:
        return n
    t = n[0]
    if t == 'c':
        return n
    if t == 'v':
        r = find(n[1])
        if r in val or r in stop:
            return ('v', r)
        d = defc.get(r)
        if d is None or len(d) != 1:
            return ('v', r)
        return expand(d[0][1], depth-1, stop)
    if t == 'neg':
        return ('neg', expand(n[1], depth, stop))
    return (t, expand(n[1], depth, stop), expand(n[2], depth, stop))

if __name__ == '__main__':
    print('classes with a unique non-copy definition: %d' %
          sum(1 for r, d in defc.items() if len(d) == 1))
    print('classes with >1 non-copy definitions: %d' %
          sum(1 for r, d in defc.items() if len(d) > 1))
    real = [t for t in lift if t[4] == Pval]
    print('lift atoms with M == P exactly: %d' % len(real))
    print()
    for k, (aid, Rw, mv, uv, M) in enumerate(real[:6]):
        print('--- lift atom %d : x%d  =  P * x%d ;  ndefs(u)=%d' % (aid, find(Rw), find(uv), ndefs(uv)))
        print('    R expanded d=3 :', pp(expand(('v', Rw), 3))[:600])
        print()
    # how free are the handles u?
    hs = Counter(ndefs(uv) for (aid, Rw, mv, uv, M) in real)
    print('handle u: #definitions histogram', dict(hs))
    hr = Counter(ndefs(Rw) for (aid, Rw, mv, uv, M) in real)
    print('R wire  : #definitions histogram', dict(hr))
