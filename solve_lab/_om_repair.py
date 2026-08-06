#!/usr/bin/env python3
"""Forward-propagate on the atom DAG, perturb the free inputs feeding the two broken
checks, and follow the repair closure: which checks break, which free inputs can absorb."""
import pickle, sys, random
from collections import defaultdict, deque
import heal_harness as H
import _om_parse as OP

D = pickle.load(open('_om_parsed2.pkl', 'rb')); astof = D['astof']
G = pickle.load(open('_om_dag.pkl', 'rb'))
usedatom = G['usedatom']; order = G['order']; checks = list(G['checks'])
avars = G['avars']; vat = G['vat']; free = set(G['free'])
checkset = set(checks)
p = H.p

# compile RHS evaluators:  x_t = RHS  from atom ('sub',('var',t),RHS)
rhs = {}
for t, k in usedatom.items():
    a = astof[k]
    assert a[0] == 'sub' and a[1] == ('var', t)
    rhs[t] = a[2]

def forward(V):
    for t in order:
        V[t] = OP.evalast(rhs[t], V)

def broken(V):
    return [k for k in checks if OP.evalast(astof[k], V) != 0]

vA = H.loadd('best_agentA_39022.json')
V0 = [0] * H.NVARS
for k, x in vA.items(): V0[k] = x
V = list(V0); forward(V)
diff = [i for i in range(H.NVARS) if V[i] != V0[i]]
print('forward() reproduces agentA exactly:', not diff, '(%d differing)' % len(diff))
B = broken(V)
print('broken checks:', len(B), B)

# --- perturbation experiment -------------------------------------------------
def try_perturb(delta_map, label):
    V = list(V0)
    for v, d in delta_map.items(): V[v] = V[v] + d
    forward(V)
    B = broken(V)
    return B

for lab, dm in [('x_7068 += 1', {7068: 1}),
                ('x_4432 += 1', {4432: 1}),
                ('x_7068,x_2964 += 1', {7068: 1, 2964: 1}),
                ('x_4432,x_24548 += 1', {4432: 1, 24548: 1}),
                ('all four += 1', {7068: 1, 2964: 1, 4432: 1, 24548: 1})]:
    B = try_perturb(dm, lab)
    print('%-24s -> %d broken: %s' % (lab, len(B), B if len(B) <= 12 else B[:12]))
