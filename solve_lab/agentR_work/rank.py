#!/usr/bin/env python3
"""Rank footprints properly.
(1) Is every cost<=6 atom a boolean-ness atom on a single variable?  If so the deliverable's
    7 is the single-atom optimum over RELATIONAL atoms and cannot be beaten with |S|=1.
(2) Are the cheap boolean atoms' variables selectors?  Relaxing a selector off {0,1} is a real
    lever (the fold could hit T with a non-boolean selection).
(3) |S|=2: can two atoms sharing equations cancel down to <=6?"""
import sys, json, re, collections, itertools
from fractions import Fraction
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
from cfgscan import E, pins, TREE

VAR = re.compile(r'x(\d+)')
atom_vars = [frozenset(int(m) for m in VAR.findall(s)) for s in E.res]
atom_eqs = collections.defaultdict(set)
for j, rr in enumerate(E.eqres):
    for c, i in rr:
        if c: atom_eqs[i].add(j)
occ = {i: len(atom_eqs[i]) for i in range(len(E.res))}

# (1)
cheap = sorted([a for a in occ if occ[a] <= 8], key=lambda a: occ[a])
print('=== (1) atoms with footprint cost <= 8 ===')
nrel = 0
for a in cheap:
    kind = 'BOOL(1 var)' if len(atom_vars[a]) == 1 else 'RELATIONAL(%d vars)' % len(atom_vars[a])
    if len(atom_vars[a]) > 1: nrel += 1
    if occ[a] <= 7: print('  cost %-3d atom %-6d %-20s %s' % (occ[a], a, kind, E.res[a][:55]))
rel = [a for a in occ if len(atom_vars[a]) > 1]
print('  min cost over RELATIONAL atoms = %d  (deliverable fails %d)' % (min(occ[a] for a in rel), 7))
print('  min cost over BOOLEAN atoms    = %d' % min(occ[a] for a in occ if len(atom_vars[a]) == 1))

# (2)
PINSET = {int(b) for b in pins}
print('=== (2) are the cheap boolean atoms selector variables? ===')
for a in cheap[:12]:
    if len(atom_vars[a]) != 1: continue
    v = next(iter(atom_vars[a]))
    print('  atom %-6d cost %-3d var x%-6d  selector/conditional-pin? %s  tree=%s'
          % (a, occ[a], v, v in PINSET, TREE.get(v)))

# (3)
print('=== (3) |S|=2 pricing over the cheapest relational atoms ===')
pool = sorted(rel, key=lambda a: occ[a])[:60]
def maxkill(S):
    rows = []
    Ss = set(S)
    for j in set().union(*[atom_eqs[a] for a in S]):
        c = {i: k for k, i in E.eqres[j] if i in Ss}
        rows.append([c.get(a, 0) for a in S])
    best = 0
    n = len(S)
    for i in range(len(rows)):
        r = rows[i]
        if all(x == 0 for x in r): continue
        # kill direction: a proportional to (r[1], -r[0]) for n=2
        vec = (r[1], -r[0]) if n == 2 else None
        if vec is None or (vec[0] == 0 and vec[1] == 0): continue
        cnt = sum(1 for q in rows if q[0] * vec[0] + q[1] * vec[1] == 0)
        best = max(best, cnt)
    return len(rows), best
res = []
for a, b in itertools.combinations(pool, 2):
    if not (atom_eqs[a] & atom_eqs[b]): continue
    t, k = maxkill([a, b])
    res.append((t - k, t, k, a, b))
res.sort()
print('  pairs sharing equations: %d ; best floors:' % len(res))
for r in res[:15]:
    print('   floor=%-4d touched=%-4d killed=%-4d atoms %d,%d  costs %d,%d'
          % (r[0], r[1], r[2], r[3], r[4], occ[r[3]], occ[r[4]]))
json.dump({'min_cost_relational': min(occ[a] for a in rel),
           'min_cost_boolean': min(occ[a] for a in occ if len(atom_vars[a]) == 1),
           'cheap_atoms': [[a, occ[a], len(atom_vars[a]), E.res[a][:60]] for a in cheap[:40]],
           'best_pairs': [list(map(int, r)) for r in res[:30]]}, open('runs/rank.json', 'w'), indent=1)
