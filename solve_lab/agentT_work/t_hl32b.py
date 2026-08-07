#!/usr/bin/env python3
"""AUDIT T35b -- test the mixed exact/divisibility fix on the SAVED |S|=32 end state, without
paying the 450 s cold rebuild.  Loads close_M32.json, calls the patched joint_pair() on the
surviving handle-less atom, and, if it clears, dumps + verifies in F's certified-faithful parse.
"""
import os, sys, json
T = '/home/user/integer_solver/solve_lab/agentT_work'
sys.path.insert(0, T)
import t_close2w as C
import t_close2wj as J
E = C.E; SL = C.SL; p = C.p; NV = C.NV
relift = C.relift; nzcount = C.nzcount

TAG = sys.argv[1] if len(sys.argv) > 1 else 'T32mix'
vv = [0]*NV
for k, val in json.load(open(os.path.join(T, 'close_M32.json'))).items():
    vv[int(k[2:])] = int(val)
relift(vv)
r = E.run(vv)
nz = [E.res[i] for i, x in enumerate(r) if x]
base = nzcount(vv)
print('loaded |S|=32 end state: %d nonzero atoms, global %d' % (len(nz), base), flush=True)
HL = [a for a in nz if a not in SL]
assert len(HL) == 1, HL
a0 = HL[0]
ok = J.joint_pair(vv, a0, base, lambda s: print(s, flush=True))
print('joint_pair ->', ok, flush=True)
relift(vv)
r = E.run(vv)
nz = [E.res[i] for i, x in enumerate(r) if x]
print('now %d nonzero atoms:' % len(nz))
for a in nz:
    print('   ', a[:100])
json.dump({'x_%d' % i: vv[i] for i in range(NV) if vv[i]},
          open(os.path.join(T, 'close_%s.json' % TAG), 'w'))
print('dumped close_%s.json' % TAG)
