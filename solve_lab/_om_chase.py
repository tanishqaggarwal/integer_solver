#!/usr/bin/env python3
"""Chase the repair closure: grow the perturbed free-input set until only the two
target checks remain broken; record the congruence each absorbing check imposes."""
import pickle, sys
from collections import defaultdict, deque
import heal_harness as H
import _om_parse as OP

D = pickle.load(open('_om_parsed2.pkl', 'rb')); astof = D['astof']
G = pickle.load(open('_om_dag.pkl', 'rb'))
usedatom = G['usedatom']; order = G['order']; checks = list(G['checks'])
avars = G['avars']; vat = G['vat']; free = set(G['free'])
checkset = set(checks)
p = H.p
rhs = {t: astof[k][2] for t, k in usedatom.items()}
TARGETS = ['((x4432-x19964)-x28730)', '((x7068-x2099)-(7376877*x642))']

def forward(V):
    for t in order: V[t] = OP.evalast(rhs[t], V)

vA = H.loadd('best_agentA_39022.json')
V0 = [0] * H.NVARS
for k, x in vA.items(): V0[k] = x

def run(S, d):
    V = list(V0)
    for v in S: V[v] += d
    forward(V)
    out = {}
    for k in checks:
        r = OP.evalast(astof[k], V)
        if r != 0: out[k] = r
    return V, out

def status(v):
    if v in free: return 'FREE'
    return 'def ' + usedatom[v]

S = {7068, 2964}
for it in range(30):
    V, br = run(S, 1)
    new = [k for k in br if k not in TARGETS]
    print('iter %d  |S|=%d  extra broken=%d  %s' % (it, len(S), len(new), new[:4]))
    if not new: break
    # examine the first new broken check; find a free var in it not yet in S
    k = new[0]
    print('   check %s  vars:' % k)
    absorbed = None
    for v in avars[k]:
        print('      x_%-6d %-46s val0=%s' % (v, status(v), str(V0[v])[:34]))
    # heuristic: add free vars of this check that are not already in S
    add = [v for v in avars[k] if v in free and v not in S]
    if not add:
        print('   no free var to add -> STUCK'); break
    S |= set(add)
    print('   adding', add)
print()
print('final S =', sorted(S))
V, br = run(S, 1)
print('broken with S shifted by 1:', list(br.keys()))
