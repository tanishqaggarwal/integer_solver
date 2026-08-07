#!/usr/bin/env python3
"""K20: descend the mux wiring from the root's two input slots down to leaf wires.
Completes the slot decode that agent F left at 47/72 stages."""
import sys, os, json, re, collections
K = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, K)
F = '/home/user/integer_solver/solve_lab/agentF_work'
sys.path.insert(0, F)
import fold as FD
import mux as MUX            # builds Engine + the Z closure at import time
from parse import node_str
from circ2 import vars_of

E = MUX.E
P = MUX.p
D = FD.points()
leafwire = {}
for l in D['leaves']:
    leafwire[l['wx']] = ('x', l['sel'])
    leafwire[l['wy']] = ('y', l['sel'])
print('leaf wires:', len(leafwire))

defrhs = MUX.defrhs
gatedpat = re.compile(r'^\(x(\d+)\*x(\d+)\)$')

memo = {}
def support(w, depth=0):
    """set of leaf selector vars reachable as sources of wire w"""
    if w in memo: return memo[w]
    if w in leafwire:
        memo[w] = {leafwire[w][1]}; return memo[w]
    memo[w] = set()
    if depth > 40: return memo[w]
    out = set()
    for z, coef in MUX.source_of(w):
        if z == 'CONST': continue
        for kind, t in MUX.mux_terms(z):
            if kind == 'gated':
                m = gatedpat.match(t)
                if m:
                    a, b = int(m.group(1)), int(m.group(2))
                    # the value wire is the non-selector one; try both
                    for u in (a, b):
                        out |= support(u, depth + 1)
            elif kind == 'free':
                out |= support(t, depth + 1)
    memo[w] = out
    return out

roots = {'A.x': 12186, 'A.y': 16742, 'B.x': 14853, 'B.y': 24908}
res = {}
for nm, w in roots.items():
    s = support(w)
    res[nm] = sorted(s)
    print(nm, 'leaf support', len(s))
json.dump(res, open(K + '/rootsupport.json', 'w'))
ch = json.load(open(K + '/chain.json'))
sel2exp = {}
for i in range(256):
    sel2exp[ch['sel'][str(i)]] = ch['exp'][str(i)]
for nm in ('A.x', 'B.x'):
    ex = sorted(sel2exp[s] for s in res[nm])
    print(nm, 'exponents:', ex)
