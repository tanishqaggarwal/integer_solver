#!/usr/bin/env python3
"""K7: cascade with boolean-selector-first seeding order."""
import sys, os, json, time, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
from cascade import Cascade, NV, P
from cascade2 import Inc

C = Cascade()
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, val in d.items(): full[int(k[2:])] = int(val)
freeset = set(C.E.free)
bools = [u for u in sorted(freeset) if full[u] in (0, 1)]
bigs = [u for u in sorted(freeset) if full[u] not in (0, 1)]
rest = [u for u in range(NV) if u not in freeset]
print('bool free', len(bools), 'big free', len(bigs), 'defined', len(rest))

I = Inc(C)
sv = {u: full[u] for u in range(NV)}
for order, lab in [(bools + bigs + rest, 'bools-first'),
                   (bools + rest + bigs, 'bools,defined,bigs'),
                   (bools, 'bools only (rest -> 0 last)')]:
    v = I.run(sv, order, verbose=False)
    bad, nz = C.score(v)
    print('%-28s derived %5d conflicts %2d nonzero-atoms %3d failing %4d score %d'
          % (lab, I.nassign, len(I.conflicts), len(nz), len(bad), 39033 - len(bad)))
    for c in I.conflicts[:10]: print('     conflict', c[0], str(c[1])[:70])
    for i in nz[:10]: print('     nz', C.atomnames[i][:80])
    json.dump({'x_%d' % i: v[i] for i in range(NV) if v[i]}, open(K + '/cand_%s.json' % lab.split()[0].replace(',', '_'), 'w'))
