#!/usr/bin/env python3
"""Circuit DAG, recognising a definition in EITHER orientation
   x_t - RHS  |  LHS - x_t  |  LHS + x_t  |  x_t + RHS ."""
import pickle
from collections import defaultdict, deque
import heal_harness as H
import _om_parse as OP

D = pickle.load(open('_om_parsed2.pkl', 'rb')); astof = D['astof']
atoms = sorted(astof)
avars = {k: sorted(OP.astvars(astof[k])) for k in atoms}
vat = defaultdict(list)
for k in atoms:
    for v in avars[k]: vat[v].append(k)

def candidates(k):
    """list of (t, other_vars) this atom could define with coefficient +-1"""
    a = astof[k]; out = []
    if a[0] in ('sub', 'add'):
        for side, oth in ((a[1], a[2]), (a[2], a[1])):
            if side[0] == 'var':
                t = side[1]
                ov = OP.astvars(oth)
                if t not in ov: out.append((t, ov))
    return out

cand = defaultdict(list)
for k in atoms:
    for t, ov in candidates(k): cand[t].append((k, ov))

solved = set(v for v in vat if v not in cand)   # never definable -> free
usedatom = {}; order = []
pend = {t: list(ks) for t, ks in cand.items()}
progress = True
while progress:
    progress = False
    for t in list(pend):
        for k, ov in pend[t]:
            if all(v in solved for v in ov):
                solved.add(t); usedatom[t] = k; order.append(t); del pend[t]
                progress = True; break
print('never-definable (hard free inputs):', len(vat) - len(cand))
print('defined:', len(order), ' still-unsolved(cyclic):', len(pend))
free = set(v for v in vat if v not in usedatom)
print('free inputs total:', len(free))
defatoms = set(usedatom.values())
checks = [k for k in atoms if k not in defatoms]
print('CHECK atoms:', len(checks))
vA = H.loadd('best_agentA_39022.json'); V = [0] * H.NVARS
for k, x in vA.items(): V[k] = x
bad = [k for k in checks if OP.evalast(astof[k], V) != 0]
print('violated checks:', len(bad), bad)
pickle.dump({'usedatom': usedatom, 'order': order, 'checks': checks,
             'free': sorted(free), 'avars': avars, 'vat': dict(vat)},
            open('_om_dag2.pkl', 'wb'))
