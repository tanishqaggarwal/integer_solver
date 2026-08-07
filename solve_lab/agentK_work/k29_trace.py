#!/usr/bin/env python3
"""K29: provenance walk-back for the B-half gap.

drive({e3,e5}) should put 8G+32G on the root's B slot (x14853,x24908) and it does not.
This records, for every variable, which atom derived it and from what, then walks the
back-cone of x14853 looking for (a) variables that were merely SEEDED rather than derived,
and (b) any variable from ABOVE the root leaking downward."""
import sys, os, json, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
from k26_drive import drive, rootpair, C, P

ch = json.load(open(K + '/chain.json'))
D = FD.points()
bypow = {}
for i_s, e in ch['exp'].items():
    bypow[e] = (int(D['leaves'][int(i_s)]['X']), int(D['leaves'][int(i_s)]['Y']))
exp2sel = {ch['exp'][str(i)]: ch['sel'][str(i)] for i in range(256)}

ON_EXP = [3, 5]
v = drive(set(exp2sel[e] for e in ON_EXP))
trace, deps = C.trace, C.deps
A, B = rootpair(v)
pred = FD.INF
for e in ON_EXP: pred = FD.add(pred, bypow[e])
print('B measured =', B)
print('B predicted=', pred)
print('match:', B == pred)

# variables that sit ABOVE the root - if any appears in the back-cone, flow is backward
ABOVE = {24468: 'target-x cmp', 18956: 'target-y cmp', 13682: 'root out mux',
         37892: 'root out mux y', 38045: 'gated out x', 10156: 'gated out y',
         22162: 'root out x', 30213: 'root out y'}

for root in (14853, 24908):
    print('\n=== back-cone of x%d ===' % root)
    print('  derived by:', C.names[trace[root]][:100] if trace.get(root) is not None else '*** SEEDED (never derived) ***')
    seen = set(); st = [root]; seeded = []; above = []
    while st:
        u = st.pop()
        if u in seen: continue
        seen.add(u)
        if u in ABOVE: above.append(u)
        if trace.get(u) is None:
            seeded.append(u); continue
        for z in deps.get(u, ()): st.append(z)
    print('  cone size:', len(seen))
    print('  variables from above the root inside the cone:',
          [(u, ABOVE[u]) for u in above] or 'NONE (flow is forward)')
    lf = set(exp2sel.values())
    sd_leaf = [u for u in seeded if u in lf]
    sd_other = [u for u in seeded if u not in lf]
    print('  seeded leaf selectors in cone:', len(sd_leaf))
    print('  seeded NON-selector variables in cone:', len(sd_other), sd_other[:20])
    for u in sd_other[:8]:
        print('     x%d  atoms:' % u)
        for i in C.var2atoms[u][:4]:
            print('        ', C.names[i][:95])
