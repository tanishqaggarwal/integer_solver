#!/usr/bin/env python3
"""Decompose S,T into sub-gates and find which depend on bits vs slack.
S = x_33469*x_29356 - x_3558^2 ; x_29356=x_29322^2 ; T = x_27713*x_29322 - x_3558*x_1326"""
import json
from agentA_harness import (p, backward_cone, load_solution, forward, gates, definer, order)

boolset = set(json.load(open('boolbits.json'))['boolvars'])
base = load_solution('best/new_instance_partial_39013.json'); forward(base)

def show(v, name):
    allv, frees = backward_cone(v)
    b = frees & boolset; s = frees - boolset
    print(f"{name:12s} x_{v}: cone {len(allv):4d} vars | free {len(frees):3d} = {len(b):3d} bool + {len(s):3d} slack | val%p={base[v]%p if base[v] else 0!r:.30}...")
    return frees, b, s

print("=== S components ===")
for v, nm in [(35389,'S'), (33469,'x_33469'), (29356,'x_29356'), (29322,'x_29322'),
              (3558,'x_3558'), (27762,'x_27762'), (14853,'x_14853'), (12186,'x_12186'),
              (24908,'x_24908'), (16742,'x_16742')]:
    show(v, nm)
print("\n=== T components ===")
for v, nm in [(6671,'T'), (27713,'x_27713'), (1326,'x_1326'), (32680,'x_32680'), (11602,'x_11602')]:
    show(v, nm)

# Which components carry the bits?
print("\n=== bit-carrying decomposition ===")
for v, nm in [(33469,'x_33469'), (29322,'x_29322'), (3558,'x_3558'), (27713,'x_27713'), (1326,'x_1326')]:
    _, b, s = backward_cone(v)[1], None, None
    frees = backward_cone(v)[1]
    bb = frees & boolset
    print(f"{nm}: {len(bb)} bool bits {'<-- CARRIES BITS' if bb else '(slack only)'}")

# show the gate rhs of S and T and key gates
gdef = {t: gates[definer[t]][1] for t in order}
print("\n=== gate rhs ===")
for v in [35389, 6671, 33469, 29356, 29322, 3558, 27713, 1326]:
    print(f"x_{v} = {gdef.get(v,'(free/pin)')[:120]}")
