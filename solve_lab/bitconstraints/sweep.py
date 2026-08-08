#!/usr/bin/env python3
"""Pass 12: local-consistency sweep.

Propagate from many different 0/1 selector assignments and count contradictions.
Two atom sets are used:
  PRIM  -- only the primitive gates (atoms with <= 8 monomials).  The large
           "combination" atoms are random Z-linear combinations of primitives,
           so a propagation that consumes them out of order can manufacture
           spurious values; excluding them keeps every derivation sound.
  ALL   -- every atom (kept for comparison).
A bit pattern that is contradiction-free under PRIM is not excluded by any
shallow (pin / gate / copy / NOT / sum) constraint.
"""
import os, json, random, sys, time
from collections import defaultdict
from closure import load, reduce_poly, atom_vars, NVARS, propagate

HERE = os.path.dirname(os.path.abspath(__file__))
D, chain = load()
atoms = D['atoms']
PRIM = [i for i, a in enumerate(atoms) if len(a) <= 8]
print(f"primitive atoms (<=8 monomials): {len(PRIM)} / {len(atoms)}")
prim_atoms = [atoms[i] for i in PRIM]


def run(bv, alist):
    val = [None] * NVARS
    for b, x in zip(chain, bv):
        val[b] = x
    det, contra = propagate(alist, val, verbose=False)
    return det, contra


rng = random.Random(4242)
cases = [('zeros', [0]*256), ('ones', [1]*256)]
for j in range(256):
    cases.append((f'w1_{j}', [1 if t == j else 0 for t in range(256)]))
for k in range(40):
    cases.append((f'rand{k}', [rng.randrange(2) for _ in range(256)]))
for k in range(40):
    bv = [0]*256
    for j in rng.sample(range(256), 2):
        bv[j] = 1
    cases.append((f'w2_{k}', bv))

t0 = time.time()
bad = []
for label, bv in cases:
    det, contra = run(bv, prim_atoms)
    if contra:
        bad.append((label, len(contra), [PRIM[i] for i in contra][:10]))
print(f"PRIM sweep: {len(cases)} bit patterns in {time.time()-t0:.0f}s; "
      f"patterns with a contradiction: {len(bad)}")
for b in bad[:20]:
    print("   ", b)

json.dump({'n_cases': len(cases), 'n_prim_atoms': len(PRIM),
           'contradicting_patterns': bad},
          open(os.path.join(HERE, 'sweep.json'), 'w'))
print("wrote sweep.json")
