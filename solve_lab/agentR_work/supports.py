#!/usr/bin/env python3
"""Compare the DELIVERABLE's defect support with the support gs2 lands in for a
single-bit configuration.  The question the ceiling in price.py cannot answer:
is the deliverable's cheap 12-equation footprint reachable at all from a
different selector configuration?"""
import sys, os, json
FW = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, FW)
from cfgscan import run_cfg, E, PIN
from fwd import NV

def footprint(v, tag):
    r = E.run(v)
    nz = [i for i, x in enumerate(r) if x]
    S = set(nz)
    eqs = [j for j, rr in enumerate(E.eqres) if any(i in S for _, i in rr)]
    bad = set(E.score(r))
    print('%-14s atoms=%-3d %s' % (tag, len(nz), nz[:12]))
    print('%-14s eqs=%-4d failing=%-4d cancelling=%d' % ('', len(eqs), len(bad), len(eqs) - len(bad)))
    return dict(tag=tag, atoms=nz, n_atoms=len(nz), eqs=eqs, n_eqs=len(eqs),
                failing=sorted(bad), n_failing=len(bad), cancelling=len(eqs) - len(bad),
                score=39033 - len(bad))

out = {}
# 1. the deliverable itself
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v = [0] * NV
for k, x in d.items():
    i = int(k[2:]) if k.startswith('x_') else int(k)
    if i < NV: v[i] = int(x)
out['deliverable'] = footprint(v, 'deliverable')
# 2. gs2 on a single bit, and on the deliverable's own configuration
for cfg in ([24601], [24601, 2081]):
    sc, nz, ok, vv = run_cfg(cfg)
    out['gs2_' + ','.join(map(str, cfg))] = footprint(vv, 'gs2 %s' % cfg)
# overlap
D = set(out['deliverable']['atoms'])
for k, f in out.items():
    if k == 'deliverable': continue
    print('%s: atom overlap with deliverable = %d ; equation overlap = %d'
          % (k, len(D & set(f['atoms'])),
             len(set(out['deliverable']['eqs']) & set(f['eqs']))))
json.dump(out, open('runs/supports.json', 'w'), indent=1)
