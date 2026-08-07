#!/usr/bin/env python3
"""THE DECIDING CURVE: floor vs gap.
Elimination degree for the backward solve is ~2^(gap between the two relaxed ladder indices).
A small gap is solvable; a large gap is not.  So: for each gap, what is the BEST floor available?
If the best floor at small gap is <= 39,026 the lever is priced out."""
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
print('selector boolean-ness atoms placed on the ladder: %d' % len(SEL))
best_by_gap = {}
for (a1, v1, c1, i1), (a2, v2, c2, i2) in itertools.combinations(SEL, 2):
    g = abs(i1 - i2)
    f = 39033 - len(atom_eqs[a1] | atom_eqs[a2])
    if g not in best_by_gap or f > best_by_gap[g][0]:
        best_by_gap[g] = (f, v1, v2, i1, i2)
print('\ngap  best-floor   pair                       elimination degree')
run = -1
for g in sorted(best_by_gap):
    f, v1, v2, i1, i2 = best_by_gap[g]
    if g <= 40 or f > 39026:
        mark = '   *** BEATS 39,026 ***' if f > 39026 else ''
        print('%-4d %-11d x%-6d(%3d)+x%-6d(%3d)   2^%-4d%s' % (g, f, v1, i1, v2, i2, g, mark))
cum = -1; curve = []
for g in sorted(best_by_gap):
    cum = max(cum, best_by_gap[g][0]); curve.append((g, cum))
print('\nBEST FLOOR REACHABLE WITH GAP <= G (i.e. with elimination degree <= 2^G):')
for G in (1, 2, 4, 8, 16, 24, 32, 48, 59, 64, 96, 128, 255):
    v = max([c for g, c in curve if g <= G], default=None)
    print('   gap <= %-4d -> best floor %s%s' % (G, v, '   BEATS' if v and v > 39026 else ''))
json.dump({'best_by_gap': {str(g): list(map(int, v)) for g, v in sorted(best_by_gap.items())}},
          open('runs/tradeoff.json', 'w'), indent=1)
