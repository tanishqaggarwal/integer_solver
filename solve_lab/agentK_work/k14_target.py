#!/usr/bin/env python3
"""K14: trace above the root gate to the target constant."""
import sys, os, json, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from cascadep import CascadeP, NV, P
from parse import node_str
from circ2 import vars_of

C = CascadeP()
defnode = {}
for a in C.E.order:
    c = C.E.cls[a]; defnode[c[1]] = c[2]
freeset = set(C.E.free)

todo = [38045, 10156, 608, 38085, 22978, 24530, 15574, 15029, 36202, 15298, 34606, 5647]
seen = set()
while todo:
    w = todo.pop(0)
    if w in seen: continue
    seen.add(w)
    n = defnode.get(w)
    print('x%d := %s' % (w, node_str(n)[:150] if n is not None else 'FREE'))
    for i in C.var2atoms[w]:
        print('        atom:', C.names[i][:150])
