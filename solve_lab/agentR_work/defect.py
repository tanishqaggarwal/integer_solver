#!/usr/bin/env python3
"""For a selector configuration: what is the DEFECT FOOTPRINT?
   - which residual atoms are nonzero after chain repair
   - how many equations those atoms occupy, and how many of those are currently satisfied
The deliverable's footprint is 7 atoms in 12 equations, 5 satisfied -> 7 failing.
A configuration whose footprint is smaller is a candidate to beat 39,026."""
import sys, os, json, time, collections
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentF_work')
import cfgscan
from cfgscan import run_cfg, E
from fwd import NV

def footprint(bits):
    sc, nz, ok, v = run_cfg(bits)
    r = E.run(v)
    nzi = [i for i, x in enumerate(r) if x]
    eqs = set()
    for j, rr in enumerate(E.eqres):
        if any(i in set(nzi) for _, i in rr): eqs.add(j)
    bad = set(E.score(r))
    return dict(bits=bits, score=sc, n_atoms=len(nzi),
                atoms=[E.res[i][:80] for i in nzi],
                n_eqs_touched=len(eqs), n_failing=len(bad),
                n_satisfied_in_footprint=len(eqs - bad))

if __name__ == '__main__':
    out = {}
    cfgs = [[24601, 2081], [24601], [2081], [47], [91], [24601, 2081, 47]]
    for c in cfgs:
        f = footprint(c)
        out[','.join(map(str, c))] = f
        print('%-18s score %d  atoms %d  eqs touched %d  failing %d  sat-in-footprint %d'
              % (str(c), f['score'], f['n_atoms'], f['n_eqs_touched'], f['n_failing'],
                 f['n_satisfied_in_footprint']), flush=True)
        json.dump(out, open('runs/defect.json', 'w'), indent=1)
