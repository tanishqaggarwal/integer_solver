#!/usr/bin/env python3
"""Pass 5: what exactly breaks at bits = all-zero?  Trace the contradiction back
to an OR-tree over the selectors and confirm it as a cardinality constraint."""
import json, os, pickle
from collections import defaultdict
from closure import load, propagate, eval_eq, atom_vars, reduce_poly, NVARS
from scan import fmt

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain = load()
BITS = set(chain)
atoms = D['atoms']
R = json.load(open(os.path.join(HERE, 'closure_results.json')))

print("=== all-zeros run ===")
z = R['zeros']
print("contradicting atoms:", z['contra_atoms'])
print("violated equations :", z['violated'])
for aid in z['contra_atoms']:
    print(f"\natom#{aid}: {fmt(atoms[aid])[:400]}")
for e in z['violated']:
    print(f"\neq#{e}: outer={D['eq_outer'][e]} natoms={len(D['eq_terms'][e])}")
    print("   atoms:", sorted({a for c, a in D['eq_terms'][e]}))
    print("   poly:", fmt(D['eq_poly'][e])[:500])

# --- back-trace the support of the contradicting atoms through the closure ---
val = [None] * NVARS
for b in chain:
    val[b] = 0
propagate(atoms, val, verbose=False)

# Which variables in the contradicting atoms, and how were they derived?
defs = defaultdict(list)   # var -> atoms that mention it
for i, a in enumerate(atoms):
    for v in atom_vars(a):
        defs[v].append(i)


def cone(seed, depth):
    """back-trace: variables reachable in `depth` levels of atom co-occurrence"""
    seen = set(seed)
    frontier = set(seed)
    for d in range(depth):
        nxt = set()
        for v in frontier:
            for i in defs[v]:
                nxt |= atom_vars(atoms[i])
        nxt -= seen
        seen |= nxt
        frontier = nxt
        print(f"   depth {d+1}: +{len(nxt)} vars, total {len(seen)}, "
              f"bits in cone {len(seen & BITS)}")
    return seen


for aid in z['contra_atoms']:
    print(f"\n--- back-cone of atom#{aid} ---")
    cone(atom_vars(atoms[aid]), 6)
