#!/usr/bin/env python3
"""THE LEVER: relax selectors off {0,1}.
A non-boolean selector b makes the mux acc' = acc + b*(S-acc) a free point on a line, so the
mux/collateral atoms can stay ZERO (they are satisfied by choosing acc' accordingly).  The only
atoms forced nonzero are the boolean-ness atoms b*(1-b) of the relaxed selectors.
Two relaxed selectors = 2 free parameters against the 2 coordinates of T -> generically solvable
for ANY boolean choice of the other 254.  Cost = |union of their boolean-atom equations|.
Beating 39,026 needs that union <= 6."""
import sys, json, itertools, collections, re
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
from cfgscan import E, pins, TREE

VAR = re.compile(r'x(\d+)')
atom_vars = [frozenset(int(m) for m in VAR.findall(s)) for s in E.res]
atom_eqs = collections.defaultdict(set)
for j, rr in enumerate(E.eqres):
    for c, i in rr:
        if c: atom_eqs[i].add(j)
PINSET = {int(b) for b in pins}

# every boolean-ness atom, tagged with its variable and whether that variable is a selector
BOOL = []
for a, vs in enumerate(atom_vars):
    if len(vs) == 1 and re.fullmatch(r'\(?x\d+\s*\*?.*', E.res[a]) and len(atom_eqs[a]):
        v = next(iter(vs))
        s = E.res[a].replace(' ', '')
        if ('*(1-x' in s) or ('*(x' in s and '-1)' in s) or ('*x' in s and ')-x' in s):
            BOOL.append((a, v, len(atom_eqs[a]), v in PINSET))
BOOL.sort(key=lambda r: r[2])
sel = [b for b in BOOL if b[3]]
print('boolean-ness atoms: %d   of which selector/conditional-pin: %d' % (len(BOOL), len(sel)))
print('cheapest selector boolean atoms:')
for a, v, c, _ in sel[:12]: print('   atom %-6d x%-6d cost %d  tree=%s' % (a, v, c, TREE.get(v)))

print('\n--- singles (1 relaxed selector: 1 parameter vs 2 conditions, needs luck) ---')
for a, v, c, _ in sel[:6]: print('   x%-6d -> %d failing -> score %d' % (v, c, 39033 - c))

print('\n--- PAIRS (2 relaxed selectors: 2 parameters vs 2 conditions, generically solvable) ---')
best = []
for (a1, v1, c1, _), (a2, v2, c2, _) in itertools.combinations(sel[:25], 2):
    u = atom_eqs[a1] | atom_eqs[a2]
    best.append((len(u), len(atom_eqs[a1] & atom_eqs[a2]), v1, v2, a1, a2, c1, c2))
best.sort()
for r in best[:12]:
    print('   x%-6d + x%-6d : union %-3d (overlap %d)  -> score %d %s'
          % (r[2], r[3], r[0], r[1], 39033 - r[0], '  *** BEATS 39026 ***' if r[0] < 7 else
             ('  (ties)' if r[0] == 7 else '')))
json.dump({'n_bool': len(BOOL), 'n_sel_bool': len(sel),
           'cheapest_sel': [[a, v, c] for a, v, c, _ in sel[:20]],
           'pairs': [list(map(int, r)) for r in best[:40]]}, open('runs/relax.json', 'w'), indent=1)
