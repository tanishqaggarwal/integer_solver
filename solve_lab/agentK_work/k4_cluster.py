#!/usr/bin/env python3
"""K4: locate the 7-atom conflict cluster. Which free inputs drive it? Are its wires stage wires?"""
import sys, os, json, collections
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from fwd import Engine, NV, compile_node
from circ2 import vars_of
from parse import node_str

E = Engine()
roles = json.load(open(F + '/stage_roles.json'))
tree = json.load(open(F + '/tree96.json'))

CLUSTER = [642, 2099, 1329, 4432, 7068, 7075, 8731, 9118, 9413, 10903, 17325, 17499,
           19964, 22665, 28599, 28730, 28961, 29854, 31864]

# stage wire membership
stagewire = {}
for g, rs in roles.items():
    r = rs[0]
    for nm in ('out', 'inA', 'inB'):
        for w in r[nm]:
            stagewire[w] = (g, nm)
print('cluster vars that are stage wires:')
for w in CLUSTER:
    print('  x%d' % w, stagewire.get(w), 'free' if w in set(E.free) else 'defined')

freeset = set(E.free)
defrhs = {E.cls[a][1]: E.atoms[a] for a in E.order}

# transitive free-input cone of each cluster var
def cone(w, cap=100000):
    seen = set(); stack = [w]; fr = set()
    while stack:
        u = stack.pop()
        if u in seen: continue
        seen.add(u)
        if u in freeset: fr.add(u); continue
        n = defrhs.get(u)
        if n is None: fr.add(u); continue
        for z in vars_of(n):
            if z != u and z not in seen: stack.append(z)
    return fr

for w in CLUSTER:
    c = cone(w)
    print('x%d cone free-inputs %d' % (w, len(c)), sorted(c)[:12])
