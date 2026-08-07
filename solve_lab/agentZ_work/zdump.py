#!/usr/bin/env python3
"""Dump the equations whose variables are all selectors, plus the multi-selector census."""
import os, sys, json, pickle, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zparse import parse, varset, reduce_L, atoms_of
from zatoms import poly, pkey

HERE = os.path.dirname(os.path.abspath(__file__))
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])

lines = [l.strip().rsplit('=', 1)[0] for l in open(EQ) if l.strip()]
IDX = [7153, 11456, 12810, 13027, 13494, 13905, 16174, 16757, 18154, 19501, 20780, 30248, 33304]

def fmt(p):
    parts = []
    for mon, c in sorted(p.items()):
        if not mon:
            parts.append("%+d" % c)
        else:
            parts.append("%+d*%s" % (c, "*".join("s%d" % v for v in mon)))
    return " ".join(parts)

for i in IDX:
    E = parse(lines[i])
    L, note = reduce_L(E)
    p = poly(L)
    print("=" * 100)
    print("EQ", i, " raw len", len(lines[i]), " vars", sorted(varset(L)))
    print("  L (expanded, s=selector) :", fmt(p))
    print("  atoms:")
    for c, a in atoms_of(L):
        print("     %+d * [%s]" % (c, fmt(poly(a))))
