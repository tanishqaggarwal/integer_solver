#!/usr/bin/env python3
"""K13: the sole surviving mod-p condition is the root chord check. Dump the root values,
find where the root output pair is pinned (the target)."""
import sys, os, json
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from cascadep import CascadeP, NV, P
from fwd import compile_node
from parse import node_str
from circ2 import vars_of

KK = 97553848499418123410591666447050222001188385549510401465815187079080512838891
C = CascadeP()
d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, val in d.items(): full[int(k[2:])] = int(val)
vc = json.load(open(K + '/varclass.json'))
order = vc['handles'] + vc['bools'] + vc['others'] + [u for u in range(NV) if u not in set(C.E.free)]
seed = {u: 0 for u in vc['handles']}
for u in vc['bools'] + vc['others']: seed[u] = full[u]
v, der = C.close(seed, order)
names = ['a_x=x12186', 'a_y=x16742', 'b_x=x14853', 'b_y=x24908', 'o_x=x22162', 'o_y=x30213']
ws = [12186, 16742, 14853, 24908, 22162, 30213]
for nm, w in zip(names, ws):
    print(nm, '=', v[w])
ax, ay, bx, by, ox, oy = [v[w] for w in ws]
print('b_x-a_x =', (bx - ax) % P, ' b_y-a_y =', (by - ay) % P)
print('root check x35389 =', v[35389])
print('recomputed:', ((ox + ax + bx + KK) * pow(bx - ax, 2, P) - pow(by - ay, 2, P)) % P)

print()
print('--- atoms pinning the root output pair ---')
for w in (22162, 30213):
    print('x%d in %d atoms:' % (w, len(C.var2atoms[w])))
    for i in C.var2atoms[w]:
        print('   ', C.names[i][:140])
print()
print('--- atoms pinning the root input wires ---')
for w in (12186, 16742, 14853, 24908):
    print('x%d in %d atoms:' % (w, len(C.var2atoms[w])))
    for i in C.var2atoms[w]:
        print('   ', C.names[i][:140])
