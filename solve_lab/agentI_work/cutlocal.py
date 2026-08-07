#!/usr/bin/env python3
"""Targeted cut search: candidates are the atoms on the defect path (the
derivation cone of a17810/a17813/a17816 and the variables they touch), ranked by
equation count.  A cut whose removal leaves the rest consistent absorbs the whole
defect at a cost equal to the number of equations it touches; below 7 beats 39,026.
"""
import os, sys, json, time, collections, itertools
from cutscan import Cutter
HERE = os.path.dirname(os.path.abspath(__file__))

C = Cutter(); M = C.M
base_val, base_conf = C.run()
base_und = sum(1 for x in base_val if x is None)
print("baseline conflicts:", base_conf, " undetermined:", base_und, flush=True)

# variables in the derivation cone of the three violated atoms
seed = set()
for a in base_conf:
    seed |= M.avars[a]
KEY = seed | {37758, 11150, 25739, 15298, 2287, 21889, 25156, 35389, 6671,
              22162, 30213, 12186, 16742, 14853, 24908, 24453}
var2atoms = collections.defaultdict(list)
for i, vs in enumerate(M.avars):
    for x in vs:
        var2atoms[x].append(i)
cands = sorted({a for v in KEY for a in var2atoms[v]},
               key=lambda a: len(M.atom_eqs[a]))
print("targeted candidates:", len(cands), flush=True)
res = {}
t0 = time.time()
out = os.path.join(HERE, 'cutlocal.json')
best = (len(base_conf), None)
# singles
for k, a in enumerate(cands):
    val, conf = C.run(disable=(a,))
    und = sum(1 for x in val if x is None)
    genuine = (len(conf) == 0 and und <= base_und + 1)
    res[str(a)] = {'nconf': len(conf), 'und': und, 'eqs': len(M.atom_eqs[a]),
                   'genuine': genuine, 'src': M.src[a][:90]}
    tag = 'GENUINE HIT' if genuine else ('partial' if len(conf) < len(base_conf) else '')
    if tag:
        print(f"  {tag} a{a} nconf={len(conf)} und={und} eqs={len(M.atom_eqs[a])} "
              f"{M.src[a][:70]}", flush=True)
    if k % 20 == 0:
        json.dump(res, open(out, 'w'))
        print(f"  {k}/{len(cands)} t={time.time()-t0:.0f}s", flush=True)
json.dump(res, open(out, 'w'))

# pairs/triples over the atoms that at least reduced the conflict count
red = [int(a) for a, r in res.items() if r['nconf'] < len(base_conf)]
print("atoms that reduce the conflict count:", len(red), red[:30], flush=True)
hits = []
for r in (2, 3):
    for combo in itertools.combinations(red, r):
        val, conf = C.run(disable=combo)
        und = sum(1 for x in val if x is None)
        if len(conf) == 0 and und <= base_und + 1:
            eqs = set()
            for a in combo:
                eqs |= {e for e, _ in M.atom_eqs[a]}
            hits.append({'set': list(combo), 'eqs': len(eqs), 'und': und})
            print(f"  GENUINE CUT {combo} touches {len(eqs)} equations", flush=True)
    json.dump({'singles': res, 'cuts': hits}, open(out, 'w'))
print("done", time.time() - t0, flush=True)
