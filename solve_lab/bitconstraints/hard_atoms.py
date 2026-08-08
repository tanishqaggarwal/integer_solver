#!/usr/bin/env python3
"""Pass 14: which atoms are *hard* constraints?

An atom A is hard iff some equation's whole left-hand side is a nonzero scalar
multiple of A (after stripping the outer constant / square).  Then A = 0 is
forced by that single equation, with no reliance on the "every atom vanishes"
modelling assumption.  We then re-verify the OR-tree derivation using hard
atoms only.
"""
import os, json, pickle
from math import gcd
from functools import reduce
from collections import defaultdict
from scan import load, support, degree, fmt

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain, BITS = load()
atoms = D['atoms']
eq_poly = D['eq_poly']


def canon(items):
    if not items:
        return ()
    g = reduce(gcd, (abs(c) for m, c in items))
    items = sorted(items)
    lead = items[0][1] // g
    sign = -1 if lead < 0 else 1
    return tuple((m, sign * c // g) for m, c in items)


key2aid = {a: i for i, a in enumerate(atoms)}
hard = defaultdict(list)         # atom id -> equations that force it alone
for i, p in enumerate(eq_poly):
    k = canon(list(p))
    aid = key2aid.get(k)
    if aid is not None:
        hard[aid].append(i)
print(f"atoms that are, by themselves, the whole content of some equation: "
      f"{len(hard)} / {len(atoms)}")
print(f"equations that reduce to a single atom: "
      f"{sum(len(v) for v in hard.values())} / {len(eq_poly)}")

OR = json.load(open(os.path.join(HERE, 'ortree.json')))
c = OR['or_constraints'][0]
chainatoms = sorted(set(c['gadget_atoms']))
pins = c['pins']
print(f"\nOR-root constraint: OR of {c['n_selectors']} selectors = 1")
print(f"  gadget atoms in its derivation: {len(chainatoms)}")
missing = [a for a in chainatoms if a not in hard]
print(f"  of these, HARD (each equal to a whole equation): "
      f"{len(chainatoms)-len(missing)}; not individually hard: {len(missing)}")
if missing:
    print("   non-hard gadget atoms:", missing[:20])
for p in pins:
    print(f"  root pin atom#{p['atom']}  hard? "
          f"{p['atom'] in hard}  forcing equations {hard.get(p['atom'], [])[:6]}")

json.dump({'hard_atoms': {str(k): v for k, v in hard.items()}},
          open(os.path.join(HERE, 'hard_atoms.json'), 'w'))
print("wrote hard_atoms.json")

# also: are the 256 booleanity atoms hard?
bool_atoms = [i for i, a in enumerate(atoms)
              if len(support(a)) == 1 and degree(a) == 2 and len(a) == 2
              and next(iter(support(a))) in BITS]
hb = [i for i in bool_atoms if i in hard]
print(f"\nbooleanity atoms for the 256 selectors: {len(bool_atoms)}, "
      f"individually hard: {len(hb)}")
