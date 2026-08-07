#!/usr/bin/env python3
"""K12: identify which stage the 3 surviving mod-p atoms belong to, and unfold them."""
import sys, os, json, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
from cascadep import CascadeP, NV, P
from parse import node_str
from circ2 import vars_of

C = CascadeP()
roles = json.load(open(F + '/stage_roles.json'))
prof = dict((g, n) for n, g in json.load(open(F + '/stage_profile.json')))
stagewire = {}
for g, rs in roles.items():
    r = rs[0]
    for nm in ('out', 'inA', 'inB'):
        for w in r[nm]: stagewire[w] = (g, nm, prof.get(g))

defnode = {}
for a in C.E.order:
    c = C.E.cls[a]; defnode[c[1]] = c[2]
freeset = set(C.E.free)

seen = set()
def unfold(w, depth=0, maxd=4):
    pre = '  ' * depth
    tag = stagewire.get(w)
    if w in freeset:
        print('%sx%d = FREE %s' % (pre, w, tag or ''))
        return
    n = defnode.get(w)
    if n is None:
        print('%sx%d = ?' % (pre, w)); return
    print('%sx%d := %s  %s' % (pre, w, node_str(n)[:110], tag or ''))
    if depth >= maxd: return
    for z in sorted(vars_of(n)):
        if z != w: unfold(z, depth + 1, maxd)

for w in [35389, 2287, 21889, 25156]:
    print('########## x%d  stagewire=%s' % (w, stagewire.get(w)))
    unfold(w)
    print()
