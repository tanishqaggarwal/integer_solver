#!/usr/bin/env python3
"""K5: print the exact local system around the 7 nonzero atoms."""
import sys, os, json, collections
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from fwd import Engine, NV, compile_node
from circ2 import vars_of
from parse import node_str

E = Engine()
freeset = set(E.free)
defrhs = {E.cls[a][1]: E.atoms[a] for a in E.order}
defnode = {}
for a in E.order:
    c = E.cls[a]
    defnode[c[1]] = c[2]

d = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
full = [0] * NV
for k, v in d.items():
    full[int(k[2:])] = int(v)

def show(w, depth=0, seen=None):
    if seen is None: seen = set()
    pre = '  ' * depth
    if w in freeset:
        print('%sx%d = FREE  (deliverable=%s)' % (pre, w, full[w]))
        return
    n = defnode.get(w)
    if n is None:
        print('%sx%d = ??' % (pre, w)); return
    s = node_str(n)
    print('%sx%d := %s   (deliverable=%s)' % (pre, w, s[:160], str(full[w])[:60]))
    if depth >= 3 or w in seen: return
    seen.add(w)
    for z in sorted(vars_of(n)):
        if z != w: show(z, depth + 1, seen)

for w in [22665, 28961, 28599, 17499, 7075, 2099, 19964]:
    print('==========', w)
    show(w)
    print()

print('=== constant values at deliverable ===')
for w in [22665, 28961, 28599, 17499, 7075, 2099, 19964, 1329, 10903, 17325, 9413, 7068, 4432, 8731, 9118, 29854, 31864, 642, 28730]:
    print('x%d = %s' % (w, full[w]))
