#!/usr/bin/env python3
"""agent AF, step 4: alias closure, locate P and every atom that carries it."""
import sys, os, pickle
from collections import Counter, defaultdict
sys.setrecursionlimit(100000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from af2_atoms import shape_of, varsof

def main():
    D = pickle.load(open(os.path.join(HERE, 'af_atoms.pkl'), 'rb'))
    E = pickle.load(open(os.path.join(HERE, 'af_defs.pkl'), 'rb'))
    atoms = D['atoms']; defs = E['defs']; consts = E['consts']

    Pval = 115792089237316195423570985008687907853269984665640564039457584007908834671663
    print('P == 2^256-2^32-977 :', Pval == 2**256 - 2**32 - 977)
    Pvar = [v for v, c in consts.items() if c == Pval][0]
    print('P is x_%d' % Pvar)

    # ---- alias closure: v defined as ('v', w) means v == w
    parent = {}
    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    ncopy = 0
    for v, lst in defs.items():
        for aid, rhs in lst:
            if rhs[0] == 'v':
                union(v, rhs[1]); ncopy += 1
    print('pure copy atoms (V - V): %d' % ncopy)
    Proot = find(Pvar)
    Pals = set(v for v in list(parent) if find(v) == Proot)
    print('P aliases (incl. itself): %d' % len(Pals))

    # constants propagated through aliases
    cval = {}
    for v, c in consts.items():
        cval[find(v)] = c
    def constof(v):
        return cval.get(find(v))

    # ---- every atom whose top level is  X - (something involving a P-alias)
    # find atoms containing any P alias variable
    hits = []
    for aid, a in enumerate(atoms):
        vs = varsof(a, set())
        if vs & Pals:
            hits.append(aid)
    print('atoms containing a P alias: %d' % len(hits))
    sh = Counter(shape_of(atoms[a]) for a in hits)
    for s, k in sh.most_common(20):
        print('   %6d  %s' % (k, s))
    pickle.dump({'Pvar': Pvar, 'Pval': Pval, 'Pals': Pals,
                 'parent': parent, 'cval': cval, 'Phits': hits},
                open(os.path.join(HERE, 'af_P.pkl'), 'wb'), 2)

if __name__ == '__main__':
    main()
