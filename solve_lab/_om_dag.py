#!/usr/bin/env python3
"""Atom-level circuit structure: which atoms DEFINE a variable, which are pure CHECKS."""
import pickle, sys
from collections import defaultdict, deque
import heal_harness as H
import _om_parse as OP

D = pickle.load(open('_om_parsed2.pkl', 'rb'))
astof = D['astof']
atoms = sorted(astof)
avars = {k: sorted(OP.astvars(astof[k])) for k in atoms}
vat = defaultdict(list)
for k in atoms:
    for v in avars[k]: vat[v].append(k)

def defines(k):
    """if atom ast is  ('sub', ('var',t), RHS)  with t not in RHS -> t"""
    a = astof[k]
    if a[0] == 'sub' and a[1][0] == 'var':
        t = a[1][1]
        if t not in OP.astvars(a[2]): return t
    return None

cand = {}
for k in atoms:
    t = defines(k)
    if t is not None: cand.setdefault(t, []).append(k)

# topological assignment: a variable is 'solved' once one of its defining atoms has all
# RHS vars solved.  Free inputs = vars with no defining atom.
allvars = set()
for k in atoms: allvars |= set(avars[k])
print('vars used in atoms:', len(allvars), 'of', H.NVARS)
free0 = [v for v in allvars if v not in cand]
print('vars with NO defining atom (free inputs):', len(free0))

solved = set(free0)
usedatom = {}
pend = deque()
remaining = dict(cand)
changed = True
order = []
while changed:
    changed = False
    for t, ks in list(remaining.items()):
        for k in ks:
            if all(v in solved or v == t for v in avars[k]):
                solved.add(t); usedatom[t] = k; order.append(t)
                del remaining[t]; changed = True; break
print('vars defined by a definable atom (topological):', len(order))
print('vars still undefined (cyclic):', len(remaining))
defatoms = set(usedatom.values())
checks = [k for k in atoms if k not in defatoms]
print('CHECK atoms (not used as a definition):', len(checks))
vA = H.loadd('best_agentA_39022.json'); V = [0] * H.NVARS
for k, x in vA.items(): V[k] = x
bad = [k for k in checks if OP.evalast(astof[k], V) != 0]
print('violated checks at agentA:', len(bad))
for k in bad: print('   ', k, '=', OP.evalast(astof[k], V))
import collections
kinds = collections.Counter()
for k in checks:
    a = astof[k]
    kinds[('vars=%d' % len(avars[k]))] += 1
print('check shapes:', kinds.most_common(10))
pickle.dump({'usedatom': usedatom, 'order': order, 'checks': checks,
             'free': free0, 'avars': avars, 'vat': dict(vat)}, open('_om_dag.pkl', 'wb'))
print('sample checks:')
for k in checks[:25]: print('   ', k)
