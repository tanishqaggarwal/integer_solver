#!/usr/bin/env python3
"""K24 END-TO-END VALIDATION: drive the real equations with EVERY selector on, close them
mod p, and compare the root's two input pairs against the group-fold prediction.
This also settles the side of leaf exponent 163."""
import sys, os, json, time
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from cascadep import CascadeP, NV, P

C = CascadeP()
vc = json.load(open(K + '/varclass.json'))
handles, bools, others = vc['handles'], vc['bools'], vc['others']
defvars = [u for u in range(NV) if u not in set(C.E.free)]
ORDER = handles + bools + others + defvars
seed = {u: 0 for u in handles}
for u in bools: seed[u] = 1
for u in others: seed[u] = 0
t0 = time.time()
v, _ = C.close(seed, ORDER)
print('closed %.1fs' % (time.time() - t0))
S = FD.SHIFT
A = ((v[12186] + S) % P, v[16742])
B = ((v[14853] + S) % P, v[24908])
print('circuit A =', A)
print('circuit B =', B)
print('A on curve:', (A[1] ** 2 - pow(A[0], 3, P) - FD.B) % P == 0)
print('B on curve:', (B[1] ** 2 - pow(B[0], 3, P) - FD.B) % P == 0)

D = FD.points()
ch = json.load(open(K + '/chain.json'))
bypow = {}
for i_s, e in ch['exp'].items():
    bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
rs = json.load(open(K + '/rootsupport.json'))
sel2exp = {ch['sel'][str(i)]: ch['exp'][str(i)] for i in range(256)}
IA = set(sel2exp[s] for s in rs['A.x']) | set(sel2exp[s] for s in rs['A.y'])
IB = set(sel2exp[s] for s in rs['B.x']) | set(sel2exp[s] for s in rs['B.y'])
MISS = 163

def foldexp(es):
    R = FD.INF
    for e in es: R = FD.add(R, bypow[e])
    return R

for lab, ia, ib in [('163->A', IA | {MISS}, IB), ('163->B', IA, IB | {MISS})]:
    fa, fb = foldexp(sorted(ia)), foldexp(sorted(ib))
    print('%-8s |IA|=%d |IB|=%d   A match: %s   B match: %s' % (lab, len(ia), len(ib), fa == A, fb == B))
