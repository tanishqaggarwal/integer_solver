#!/usr/bin/env python3
"""Minimum-cost cut.

Established: freeing X24453 (the K pin, atom a40368, ONE equation) lets A=0 be
solved, and freeing any one coordinate in B lets B=0 be solved; the mod-p system
then closes completely (0 conflicts).  Cost = number of equations whose atom
combination is nonzero.  a40368 costs 1; the question is the cheapest second cut.

This script enumerates every atom in the derivation chain of x1,y1,x2,y2,x3,y3,
ranked by equation count, and tests each as the second cut.
"""
import os, json, time, collections
from cutscan import Cutter
from fp import P
HERE = os.path.dirname(os.path.abspath(__file__))

C = Cutter(); M = C.M
val, conf = C.run()
base_und = sum(1 for x in val if x is None)
COORD = {'x1': 12186, 'y1': 16742, 'x2': 14853, 'y2': 24908,
         'x3': 22162, 'y3': 30213}
x1, y1, x2, y2, x3, y3, K = (val[v] for v in
                             (12186, 16742, 14853, 24908, 22162, 30213, 24453))
Kp = ((y2 - y1)**2 * pow((x2 - x1)**2, -1, P) - x3 - x1 - x2) % P
print("baseline conflicts", conf, "undetermined", base_und, flush=True)

# derivation chain atoms of the coordinates (reasons recorded by Cutter.run?)
# rebuild with boolscore to get reasons
from boolscore import Fast
F = Fast()
polwit = lambda u, r: (F.witp[u] if F.witp[u] in r else r[0])
F.run(polwit)
R = F.reason
chain_atoms = set()
for nm, v in COORD.items():
    seen = set(); st = [v]
    while st:
        u = st.pop()
        if u in seen:
            continue
        seen.add(u)
        r = R[u]
        if r in (None, 'dec', 'pre'):
            continue
        chain_atoms.add(r)
        st.extend(y for y in M.avars[r] if y != u)
cands = sorted(chain_atoms, key=lambda a: len(M.atom_eqs[a]))
print("chain atoms:", len(cands),
      "cheapest:", [(a, len(M.atom_eqs[a])) for a in cands[:8]], flush=True)

results = []
t0 = time.time()
for k, a in enumerate(cands):
    cost = 1 + len(M.atom_eqs[a])
    if cost >= 7:
        break                      # cannot beat the deliverable
    # what does releasing a free?  find which coordinate becomes undetermined
    v1, c1 = C.run(disable=(40368, a), preassign={24453: Kp})
    freed = [nm for nm, vv in COORD.items() if v1[vv] is None]
    results.append({'atom': a, 'eqs': len(M.atom_eqs[a]), 'cost_bound': cost,
                    'freed': freed, 'nconf': len(c1),
                    'src': M.src[a][:90]})
    print(f"  a{a} eqs={len(M.atom_eqs[a])} costbound={cost} freed={freed} "
          f"nconf={len(c1)}  {M.src[a][:70]}", flush=True)
json.dump(results, open(os.path.join(HERE, 'mincut.json'), 'w'), indent=1)
print("candidates with cost bound < 7:", len(results), " t=%.0f" % (time.time() - t0))
