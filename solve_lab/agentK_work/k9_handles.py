#!/usr/bin/env python3
"""K9: detect handle variables (all their atom coefficients divisible by p), seed them FIRST at 0,
then booleans, then let the cascade derive everything else."""
import sys, os, json, collections, time
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
from cascade import Cascade, NV, P
from cascade2 import Inc

C = Cascade()
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, val in d.items(): full[int(k[2:])] = int(val)
freeset = set(C.E.free)

# structural handle detection: coefficient of u in every atom containing u is divisible by p.
# Evaluate coefficients at the deliverable assignment (coefficients here are products of
# constants/other vars, so use the deliverable point; then re-check at a second point).
import random
def coeffs_at(base):
    out = {}
    for u in sorted(freeset):
        cs = []
        for i in C.var2atoms[u]:
            old = base[u]
            base[u] = old + 1; c1 = C.evala(i, base)
            base[u] = old; c0 = C.evala(i, base)
            cs.append(c1 - c0)
        out[u] = cs
    return out

cf = coeffs_at(list(full))
handles = [u for u, cs in cf.items() if cs and all(c != 0 and c % P == 0 for c in cs)]
print('handle-like free vars:', len(handles))
nonlin = [u for u, cs in cf.items() if not cs]
print('free vars in no atom:', len(nonlin))
bools = [u for u in sorted(freeset) if u not in set(handles) and full[u] in (0, 1)]
others = [u for u in sorted(freeset) if u not in set(handles) and full[u] not in (0, 1)]
print('bools', len(bools), 'other free', len(others))
json.dump({'handles': handles, 'bools': bools, 'others': others}, open(K + '/varclass.json', 'w'))

I = Inc(C)
sv = {u: full[u] for u in range(NV)}
sv0 = dict(sv)
for u in handles: sv0[u] = 0

for order, seedv, lab in [
    (handles + bools + others + [u for u in range(NV) if u not in freeset], sv0, 'handles0,bools,others'),
    (handles + bools + others + [u for u in range(NV) if u not in freeset], sv, 'handlesDeliv,bools,others'),
]:
    v = I.run(seedv, order, verbose=False)
    bad, nz = C.score(v)
    print('%-30s derived %5d conflicts %2d nonzero-atoms %3d failing %4d score %d'
          % (lab, I.nassign, len(I.conflicts), len(nz), len(bad), 39033 - len(bad)))
    for c in I.conflicts[:12]: print('     conflict', c[0], str(c[1])[:70])
    for i in nz[:12]: print('     nz', C.atomnames[i][:80])
    json.dump({'x_%d' % i: v[i] for i in range(NV) if v[i]}, open(K + '/cand_k9_%s.json' % lab.split(',')[0], 'w'))
