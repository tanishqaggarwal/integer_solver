#!/usr/bin/env python3
"""K38: if a gadget's output is vacuous whenever its two inputs COINCIDE, then the cheapest
coincidence in the whole circuit is not two equal points at all -- it is TWO DEAD INPUTS.

All leaf selectors off => every leaf wire is pinned to 0 by its off-pin => every slot wire is
0 => at every gadget a == b == (0,0).  If a gadget's liveness gate can be held at 1 while its
inputs are dead, its residuals are vacuous by Q's mechanism and its output is free -- and a
free output at any gadget propagates up as a pass-through and can be set to the target.

Whether that is reachable depends on whether the liveness bits are FORCED to 0 by the
selectors or can be set independently.  Measure it rather than argue it: score the mod-p
system for the all-off configuration with the non-leaf booleans held 0, held 1, and derived."""
import sys, os, json, time
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from cascadep import CascadeP, NV, P
from k26_drive import FORBID

C = CascadeP()
vc = json.load(open(K + '/varclass2.json'))
handles, leafsel, otherbools, wires = vc['handles'], vc['leafsel'], vc['otherbools'], vc['wires']
defvars = [u for u in range(NV) if u not in set(C.E.free)]
ch = json.load(open(K + '/chain.json'))
exp2sel = {ch['exp'][str(i)]: ch['sel'][str(i)] for i in range(256)}
S = FD.SHIFT
TGT = (int(FD.points()['target']['X']), int(FD.points()['target']['Y']))


def run(on, ob, label, forbid=FORBID):
    seed = {u: 0 for u in handles}
    for u in leafsel: seed[u] = 1 if u in on else 0
    if ob is None:
        order = handles + leafsel + defvars + wires + otherbools
    else:
        for u in otherbools: seed[u] = ob
        order = handles + leafsel + otherbools + defvars + wires
    t0 = time.time()
    v, _ = C.close(seed, order, forbid=forbid)
    bad, nz = C.score(v)
    print('%-42s failing eqs mod p = %-6d nonzero atoms = %-5d  (%.0fs)'
          % (label, bad, len(nz), time.time() - t0))
    print('      root gate x15298 = %d   root out = (%d..., %d...)   out==target: %s'
          % (v[15298], (v[22162] + S) % P % 10 ** 9, v[30213] % 10 ** 9,
             ((v[22162] + S) % P, v[30213]) == TGT))
    return v, bad, nz


print('BASELINE for scale: the shared deliverable scores 7 failing equations over Z.')
print()
run(set(), 0, 'all selectors OFF, other bools = 0')
run(set(), 1, 'all selectors OFF, other bools = 1')
run(set(), None, 'all selectors OFF, other bools derived')
run({exp2sel[0]}, 0, 'one leaf ON (e0), other bools = 0')
run({exp2sel[e] for e in (0, 1, 3)}, 0, 'three leaves ON, other bools = 0')
