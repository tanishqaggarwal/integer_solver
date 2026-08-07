#!/usr/bin/env python3
"""K6: forward-eval repair. Free inputs from the deliverable, then try zeroing the
gate-21279 output pair (x8731,x9118) and re-solving the four handle pins."""
import sys, os, json, time
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from fwd import Engine, NV

p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
E = Engine()
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, v in d.items():
    full[int(k[2:])] = int(v)

def trial(mods, label):
    v = [0] * NV
    for i in E.free: v[i] = full[i]
    for k, val in mods.items(): v[k] = val
    r = E.run(v)
    nz = [(E.res[i], r[i]) for i, x in enumerate(r) if x]
    bad = E.score(r)
    print('%-40s nonzero-res %d  failing %d  score %d' % (label, len(nz), len(bad), 39033 - len(bad)))
    for a, val in nz[:12]:
        print('     ', a[:90], '  val/p =', val // p if val % p == 0 else 'NOT div by p (%d digits)' % len(str(abs(val))))
    return v, r, nz, bad

trial({}, 'baseline free inputs')
# step 1: zero the gate output pair and its two handles
v, r, nz, bad = trial({8731: 0, 9118: 0, 1329: 0, 10903: 0}, 'zero out-pair + 2 handles')
# step 2: also fix the two congruence pins by setting the free lhs equal to the rhs
x2099 = v[2099]; x19964 = v[19964]
v2, r2, nz2, bad2 = trial({8731: 0, 9118: 0, 1329: 0, 10903: 0,
                           17325: 0, 7068: x2099, 9413: 0, 4432: x19964}, 'all four pins repaired')
json.dump({'x_%d' % i: v2[i] for i in range(NV) if v2[i]},
          open('/home/user/integer_solver/solve_lab/agentK_work/cand_k6.json', 'w'))
print('wrote cand_k6.json')
