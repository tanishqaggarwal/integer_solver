#!/usr/bin/env python3
"""PRICE THE OBSTRUCTION before attempting any solve.

Backward-solve structure for two relaxed selectors at ladder indices i < j:
  * stages j+1..255 are boolean and the fold law is invertible -> walk the root back to a
    CONCRETE required point R_j.  Cheap, no degree growth.
  * stages 0..i-1 are boolean -> forward to a CONCRETE point A.  Cheap.
  * at i:  acc_i = A + t1*(S_i - A)   -- a LINE, degree 1 in t1.
  * stages i+1..j-1 push that line through the chord law: each stage is a quadratic-over-quadratic
    rational map, so the degree in t1 roughly DOUBLES per stage.
  * at j:  eliminate t2 by requiring R_j to be collinear with B(t1) and S_j(B(t1)) -> ONE
    polynomial equation in t1 of degree ~2^(j-i).
So the cost is set by the GAP j-i, not by depth from the root.  Exhaustive over all selector
boolean-ness atoms, not just the cheapest few."""
import sys, json, re, collections, itertools
sys.path.insert(0, '.')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
import model
from cfgscan import E, pins as PINS
P, S = model.P, model.S
lad = json.load(open('ladder.json'))['ladder']
LX = {int(e[1]): i for i, e in enumerate(lad)}
pj = json.load(open('/home/user/integer_solver/solve_lab/agentF_work/pins.json'))
def lidx(v):
    e = pj.get(str(v))
    if not e or len(e) < 2: return None
    return LX.get((int(e[1][1]) + S) % P)

VAR = re.compile(r'x(\d+)')
atom_vars = [frozenset(int(m) for m in VAR.findall(s)) for s in E.res]
atom_eqs = collections.defaultdict(set)
for j, rr in enumerate(E.eqres):
    for c, i in rr:
        if c: atom_eqs[i].add(j)
PINSET = {int(b) for b in PINS}
SEL = []
for a, vs in enumerate(atom_vars):
    if len(vs) != 1: continue
    v = next(iter(vs))
    if v not in PINSET: continue
    s = E.res[a].replace(' ', '')
    if ('*(1-x' in s) or ('*(x' in s and '-1)' in s) or ('*x' in s and ')-x' in s):
        li = lidx(v)
        if li is not None: SEL.append((a, v, len(atom_eqs[a]), li))
print('selector boolean-ness atoms with a ladder position: %d' % len(SEL))
print('cost range %d..%d' % (min(s[2] for s in SEL), max(s[2] for s in SEL)))

print('\n=== EXHAUSTIVE pair scan over ALL %d, ranked by floor then by gap ===' % len(SEL))
rows = []
for (a1, v1, c1, i1), (a2, v2, c2, i2) in itertools.combinations(SEL, 2):
    u = len(atom_eqs[a1] | atom_eqs[a2])
    if u > 8: continue
    rows.append((u, abs(i1 - i2), v1, v2, i1, i2))
rows.sort()
print('%-7s %-6s %-38s %s' % ('floor', 'gap', 'pair (ladder indices)', 'elimination degree ~2^gap'))
for u, g, v1, v2, i1, i2 in rows[:25]:
    print('%-7d %-6d x%-6d(%3d) + x%-6d(%3d)%14s 2^%d' % (39033 - u, g, v1, i1, v2, i2, '', g))
beat = [r for r in rows if 39033 - r[0] > 39026]
print('\npairs that BEAT 39,026: %d ; their minimum gap = %s'
      % (len(beat), min((r[1] for r in beat), default=None)))
ties = [r for r in rows if 39033 - r[0] == 39026]
print('pairs that TIE  39,026: %d ; their minimum gap = %s'
      % (len(ties), min((r[1] for r in ties), default=None)))
json.dump({'n_sel': len(SEL),
           'rows': [[int(39033 - u), int(g), int(v1), int(v2), int(i1), int(i2)] for u, g, v1, v2, i1, i2 in rows[:80]],
           'min_gap_beating': min((r[1] for r in beat), default=None),
           'min_gap_tying': min((r[1] for r in ties), default=None)}, open('runs/obstruct.json', 'w'), indent=1)
