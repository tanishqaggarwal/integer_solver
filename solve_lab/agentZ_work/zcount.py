#!/usr/bin/env python3
"""Agent Z: is there a COUNTING network anywhere?  A cardinality constraint on
|S| needs boolean values to be ADDED somewhere.  Census over every atom in the
instance of the shapes  y = a+b  (adder)  vs  y = a+b-ab  (OR, cannot count)."""
import os, sys, json, pickle, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zparse import parse, varset, reduce_L, atoms_of
from zatoms import poly

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])
eqp, boolvars = pickle.load(open(os.path.join(HERE, 'zbool.pkl'), 'rb'))
print("boolean vars:", len(boolvars))

sys.setrecursionlimit(100000)
lines = [l.strip().rsplit('=', 1)[0] for l in open(EQ) if l.strip()]

seen = set()
adders = []
ors = []
allbool_affine = collections.Counter()
for i, lhs in enumerate(lines):
    E = parse(lhs)
    L, _ = reduce_L(E)
    for c, a in atoms_of(L):
        p = poly(a)
        vs = set()
        for m in p:
            vs |= set(m)
        if not vs or not (vs <= boolvars):
            continue
        k = tuple(sorted(p.items()))
        if k in seen:
            continue
        seen.add(k)
        deg = max(len(m) for m in p)
        if deg == 1 and len(vs) >= 3:
            allbool_affine[len(vs)] += 1
            cs = sorted(p[(v,)] for v in vs)
            # adder shape: coefficients {+1,+1,-1} (y = a+b) up to global sign
            if len(vs) == 3 and p.get((), 0) == 0:
                vals = sorted(p[(v,)] for v in vs)
                if vals in ([-1, 1, 1], [-1, -1, 1]):
                    adders.append((i, p))
        if deg == 2:
            # OR:  y - a - b + a*b
            quad = [m for m in p if len(m) == 2]
            if len(quad) == 1 and len(vs) == 3:
                (qa, qb), = quad
                if qa != qb:
                    ors.append((i, p))

print("distinct all-boolean atoms:", len(seen))
print("all-boolean AFFINE atoms by #vars:", sorted(allbool_affine.items()))
print("ADDER-shaped atoms  (y = a+b over booleans):", len(adders))
for i, p in adders[:10]:
    print("   eq%-6d %s" % (i, sorted(p.items())))
print("OR-shaped atoms (y = a+b-ab over booleans):", len(ors))
for i, p in ors[:6]:
    print("   eq%-6d %s" % (i, sorted(p.items())))
