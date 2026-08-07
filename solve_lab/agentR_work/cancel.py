#!/usr/bin/env python3
"""Can a single-bit configuration's defect be given CANCELLATION?
The deliverable's defect sits in ~13 equations of which 6 cancel; gs2's single-bit defect sits
in 20 equations of which 0 cancel.  For each of those 20 equations, price the cancellation:
how many OTHER atoms does it contain (a live partner is needed for the sum to vanish), and what
does turning such a partner on cost elsewhere?"""
import sys, json, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
from cfgscan import run_cfg, E

sc, nz, ok, v = run_cfg([24601])
r = E.run(v)
S = {i for i, x in enumerate(r) if x}
print('live atoms', sorted(S), 'score(E)', sc, flush=True)

atom_eqs = collections.defaultdict(set)
for j, rr in enumerate(E.eqres):
    for _, i in rr: atom_eqs[i].add(j)

foot = [j for j, rr in enumerate(E.eqres) if any(i in S for _, i in rr)]
rows = []
for j in foot:
    ats = [i for _, i in E.eqres[j]]
    live = [i for i in ats if i in S]
    dead = [i for i in ats if i not in S]
    cost = sorted({i: len(atom_eqs[i]) for i in dead}.items(), key=lambda kv: kv[1])
    rows.append((j, len(ats), len(live), len(dead), cost[:3]))

solo = [x for x in rows if x[1] == 1]
print('equations in footprint: %d' % len(rows))
print('  with only ONE atom total (can NEVER cancel): %d' % len(solo))
print('  with a dead partner available              : %d' % sum(1 for x in rows if x[3] > 0))
cheap = [x for x in rows if x[3] > 0 and x[4] and x[4][0][1] <= 3]
print('  whose cheapest partner touches <=3 equations: %d' % len(cheap))
for x in rows[:20]:
    print('   eq %-6d atoms=%-3d live=%d dead=%d cheapest partners(eqs touched)=%s'
          % (x[0], x[1], x[2], x[3], x[4]))
json.dump({'footprint': [list(map(str, x[:4])) for x in rows],
           'n_eqs': len(rows), 'n_uncancellable': len(solo),
           'n_with_partner': sum(1 for x in rows if x[3] > 0),
           'n_cheap_partner': len(cheap)}, open('runs/cancel.json', 'w'), indent=1)
