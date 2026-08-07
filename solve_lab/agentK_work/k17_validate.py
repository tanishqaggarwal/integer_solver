#!/usr/bin/env python3
"""K17 VALIDATION: drive the real equations with a chosen leaf ON-set, close them mod p,
read the root input pairs out of the circuit, and compare with the group-fold prediction."""
import sys, os, json, time, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from cascadep import CascadeP, NV, P
import fold as FD

C = CascadeP()
D = FD.points()
idx = {l['sel']: (int(l['X']), int(l['Y'])) for l in D['leaves']}
SEL = sorted(idx)
vc = json.load(open(K + '/varclass.json'))
handles, bools, others = vc['handles'], vc['bools'], vc['others']
defvars = [u for u in range(NV) if u not in set(C.E.free)]
ORDER = handles + bools + others + defvars
SELSET = set(SEL)

def run(on):
    seed = {u: 0 for u in handles}
    for u in bools: seed[u] = 1 if u in on else 0
    for u in others: seed[u] = 0
    v, _ = C.close(seed, ORDER)
    return v

ROOT = dict(ax=12186, ay=16742, bx=14853, by=24908, ox=22162, oy=30213, fx=13682, fy=37892)
S = FD.SHIFT
print('%-8s %-12s %-12s %-12s %s' % ('leaf', 'A-side?', 'B-side?', 'final==leaf?', 'time'))
res = {}
t0 = time.time()
for s in SEL[:6]:
    v = run({s})
    X, Y = idx[s]
    ax, ay = (v[ROOT['ax']] + S) % P, v[ROOT['ay']]
    bx, by = (v[ROOT['bx']] + S) % P, v[ROOT['by']]
    fx, fy = (v[ROOT['fx']] + S) % P, v[ROOT['fy']]
    res[s] = ('A' if (ax, ay) == (X, Y) else '', 'B' if (bx, by) == (X, Y) else '')
    print('%-8d %-12s %-12s %-12s %.1fs' % (s, (ax, ay) == (X, Y), (bx, by) == (X, Y),
                                            (fx, fy) == (X, Y), time.time() - t0))
    if s == SEL[0]:
        print('   leaf   X,Y =', X, Y)
        print('   circuit A =', ax, ay)
        print('   circuit B =', bx, by)
        print('   circuit final =', fx, fy)
