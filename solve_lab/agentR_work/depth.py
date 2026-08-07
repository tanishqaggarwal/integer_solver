#!/usr/bin/env python3
"""Depth profile: where do the cheap relaxable selectors sit on the ladder?
Stages downstream of ladder index i = 255 - i.  The backward solve pushes
acc' = acc + t*(S-acc) symbolically through every downstream stage, so the composite degree
grows with that count: near the root (high index) small, near the leaves (low index) hopeless.

Lookup fixed: pins.json stores [[wire_y, val_y], [wire_x, val_x]] -- the SECOND entry is the
x-coordinate, in the unshifted frame, so ladder x = (val + S) mod P.  Validated below against the
four-agent-confirmed ON-set {24601 -> 72, 2081 -> 235}."""
import json, sys, collections
sys.path.insert(0, '.')
import model
P, S = model.P, model.S
lad = json.load(open('ladder.json'))['ladder']
LX = {int(e[1]): i for i, e in enumerate(lad)}
pins = json.load(open('/home/user/integer_solver/solve_lab/agentF_work/pins.json'))

def ladder_index(v):
    e = pins.get(str(v))
    if not e or len(e) < 2: return None
    return LX.get((int(e[1][1]) + S) % P)

print('VALIDATION (four agents independently confirm these):')
for v, want in ((24601, 72), (2081, 235)):
    got = ladder_index(v)
    print('   x%-6d -> ladder %s   expected %s   %s' % (v, got, want, 'OK' if got == want else 'MISMATCH'))

CAND = [(8509, 33095, 3), (7889, 19326, 6), (8510, 28825, 6), (8511, 4362, 7)]
print('\nDEPTH PROFILE of the cheap relaxable selectors:')
prof = {}
for atom, var, cost in CAND:
    i = ladder_index(var)
    prof[var] = i
    print('   atom %-5d x%-6d cost %-2d -> ladder index %-5s  stages downstream = %s'
          % (atom, var, cost, i, (255 - i) if i is not None else '?'))

print('\nPAIRS, floor vs worst (deepest) position — the solve is driven by the LOWEST index:')
import itertools, re
atom_eqs = collections.defaultdict(set)
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
from cfgscan import E
for j, rr in enumerate(E.eqres):
    for c, i in rr:
        if c: atom_eqs[i].add(j)
rows = []
for (a1, v1, c1), (a2, v2, c2) in itertools.combinations(CAND, 2):
    u = len(atom_eqs[a1] | atom_eqs[a2])
    i1, i2 = prof[v1], prof[v2]
    lo = min(x for x in (i1, i2) if x is not None) if (i1 is not None and i2 is not None) else None
    rows.append((u, lo, v1, v2, i1, i2))
rows.sort(key=lambda r: (r[0], -(r[1] if r[1] is not None else -1)))
for u, lo, v1, v2, i1, i2 in rows:
    print('   x%-6d(idx %-4s) + x%-6d(idx %-4s) : floor %d  -> deepest stage count = %s'
          % (v1, i1, v2, i2, 39033 - u, (255 - lo) if lo is not None else '?'))
json.dump({'profile': {str(k): v for k, v in prof.items()},
           'pairs': [[int(u), (int(lo) if lo is not None else -1), int(v1), int(v2)] for u, lo, v1, v2, _, _ in rows]},
          open('runs/depth.json', 'w'), indent=1)
