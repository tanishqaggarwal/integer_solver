#!/usr/bin/env python3
"""Iterate variable/equation closure and test tangent consistency at each expansion.
Does the coupling close into a finite consistent subsystem, or stay obstructed?"""
import heal_harness as H
from jac_lib import D
import flint
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
H.val[17325]=0; H.val[9413]=0; H.forward()
H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]; H.forward()
F1=set(H.fails())
print('failing',len(F1))

# eq -> free support map
def eqfree(i):
    s=set()
    for v in H.eqvars[i]: s|=H.anc.get(v,{v})
    return s & H.freeinp
# free input -> eqs that depend on it (via ancestor)
from collections import defaultdict
var2eq=defaultdict(set)
for i in range(len(H.eqcode)):
    for j in eqfree(i): var2eq[j].add(i)

cols=set()
for i in F1: cols|=eqfree(i)
for it in range(6):
    # affected eqs = all eqs touching any col
    aff=set()
    for j in cols: aff|=var2eq[j]
    newcols=set()
    for i in aff: newcols|=eqfree(i)
    print(f'iter {it}: cols={len(cols)} affected_eqs={len(aff)} -> newcols={len(newcols)}')
    if newcols==cols:
        print('CLOSED'); break
    cols=newcols
    if len(cols)>4000:
        print('coupling does NOT close (>4000 free vars) — global'); break
