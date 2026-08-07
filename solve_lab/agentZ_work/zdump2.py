#!/usr/bin/env python3
"""Dump, for a few selectors, EVERY equation that mentions them, atom by atom."""
import os, sys, json, pickle, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zparse import parse, varset, reduce_L, atoms_of
from zatoms import poly, pkey

HERE = os.path.dirname(os.path.abspath(__file__))
sel = set(json.load(open(os.path.join(HERE, 'zsel.json')))['selectors'])
D = pickle.load(open(os.path.join(HERE, 'zclass.pkl'), 'rb'))
atom_seen = D['atom_seen']

def fmt(p):
    out = []
    for mon, c in sorted(p.items(), key=lambda kv: (len(kv[0]), kv[0])):
        if not mon:
            s = "%+d" % c
        else:
            nm = "*".join(("s%d" if v in sel else "x%d") % v for v in mon)
            s = "%+d*%s" % (c, nm) if abs(c) != 1 else ("%s%s" % ('+' if c > 0 else '-', nm))
        out.append(s)
    return " ".join(out) if out else "0"

# 1. class exemplars
print("### exemplars per class")
byclass = collections.defaultdict(list)
for k, (p, tag, sv) in atom_seen.items():
    byclass[tag].append(p)
for tag, ps in sorted(byclass.items()):
    print("--", tag, len(ps))
    for p in ps[:4]:
        print("     ", fmt(p))

# 2. full equations mentioning selector s
EQ = os.path.join(HERE, '..', '..', 'EQUATIONS.txt')
lines = [l.strip().rsplit('=', 1)[0] for l in open(EQ) if l.strip()]
eq_sel = D['eq_sel_atoms']
target = sorted(sel)[0]
print()
print("### every equation mentioning selector s%d" % target)
n = 0
for i, ats in sorted(eq_sel.items()):
    svs = set()
    for c, k, tag, sv in ats:
        svs |= set(sv)
    if target not in svs:
        continue
    E = parse(lines[i]); L, _ = reduce_L(E)
    allv = varset(L)
    print("EQ %d  nvars=%d  selectors=%s" % (i, len(allv), sorted(svs)))
    for c, a in atoms_of(L):
        p = poly(a)
        sv = set()
        for m in p:
            sv |= set(m) & sel
        mark = "  <<<" if target in sv else ""
        tag = atom_seen.get(pkey(p), (None, 'NOSEL', None))[1]
        print("    %+6d * [%s]   %s%s" % (c, fmt(p), tag, mark))
    n += 1
    if n >= 8:
        break
