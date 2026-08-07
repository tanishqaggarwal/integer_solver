#!/usr/bin/env python3
"""The cheap footprints are BOOLEAN-NESS atoms x*(1-x) on a single variable.
cost(S) is only a floor: making atom a nonzero means driving its variable off {0,1},
and that variable feeds OTHER atoms.  Price the true footprint = every atom containing
that variable, and the equations they touch."""
import sys, json, re, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
from cfgscan import E

VAR = re.compile(r'x(\d+)')
atom_vars = [frozenset(int(m) for m in VAR.findall(s)) for s in E.res]
var_atoms = collections.defaultdict(set)
for a, vs in enumerate(atom_vars):
    for v in vs: var_atoms[v].add(a)
atom_eqs = collections.defaultdict(set)
for j, rr in enumerate(E.eqres):
    for c, i in rr:
        if c: atom_eqs[i].add(j)

occ = {i: len(atom_eqs[i]) for i in range(len(E.res))}
cands = sorted(occ.items(), key=lambda kv: kv[1])[:60]
rows = []
for a, c in cands:
    vs = atom_vars[a]
    if len(vs) != 1:
        rows.append((a, c, len(vs), None, None, E.res[a][:60])); continue
    v = next(iter(vs))
    coll = var_atoms[v] - {a}                     # atoms that also move when x_v leaves {0,1}
    eqs = set(atom_eqs[a])
    for b in coll: eqs |= atom_eqs[b]
    rows.append((a, c, 1, len(coll), len(eqs), E.res[a][:60]))

rows.sort(key=lambda r: (r[4] if r[4] is not None else 10 ** 9))
print('%-7s %-6s %-5s %-9s %-11s %s' % ('atom', 'cost', 'nvars', 'collateral', 'true-touch', 'expression'))
for r in rows[:30]:
    print('%-7d %-6d %-5d %-9s %-11s %s' % (r[0], r[1], r[2], r[3], r[4], r[5]), flush=True)
json.dump([{'atom': r[0], 'solo_cost': r[1], 'nvars': r[2], 'collateral_atoms': r[3],
            'true_eqs_touched': r[4], 'expr': r[5]} for r in rows],
          open('runs/reach2.json', 'w'), indent=1)
