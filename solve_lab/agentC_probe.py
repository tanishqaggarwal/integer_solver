#!/usr/bin/env python3
"""For each deep control (a=x_12186,b=x_14853,c=x_16742) and the x_24908 cone, perturb by 1, find
which currently-satisfied eqs break (root != 0 mod p) and what OTHER free inputs those eqs contain
(compensator candidates). Report the coupling structure and compensator-set sizes."""
import json, time
from agentC_common import (p, gates, order, definer, gcode, forward, partial_forward, downstream_ks,
                           val, freeinp, ns, lines, eqcode, eqvars, load_best, CORE, posof,
                           NVARS, pinned, rootcode_of)
from collections import defaultdict

gate_defs = {t: (rhs, vids) for t, rhs, vids in gates}
def freecone(root):
    seen=set(); leaves=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x][1]: st.append(u)
        elif x in freeinp: leaves.add(x)
    return leaves

best=load_best(); forward()
F0=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
def broken_by(h):
    ks=downstream_ks(h)
    aff=set(eqbyvar.get(h,()))
    for k in ks: aff|=eqbyvar.get(order[k],set())
    aff-=F0
    o=val[h]; val[h]=o+1; partial_forward(ks)
    br=[]
    for i in aff:
        if eval(rootcode_of(i),ns)%p!=0: br.append(i)
    val[h]=o; partial_forward(ks)
    return br

controls={'a=x_12186':12186,'b=x_14853':14853,'c=x_16742':16742}
xcone=sorted(freecone(24908))
print(f"x_24908 cone: {len(xcone)} free inputs")
allcomp=set()
for name,h in controls.items():
    br=broken_by(h)
    comp=set()
    for i in br: comp|=(eqvars[i]&freeinp)
    comp-={h}
    allcomp|=comp
    print(f"{name}: breaks {len(br)} sat-eqs; compensator free inputs available: {len(comp)}")
# x_24908 cone collectively
cone_comp=set()
cone_br=set()
for h in xcone:
    br=broken_by(h)
    cone_br|=set(br)
    for i in br: cone_comp|=(eqvars[i]&freeinp)
cone_comp-=set(xcone)
print(f"x_24908 cone (all): breaks {len(cone_br)} sat-eqs; compensators: {len(cone_comp)}")
allcomp|=cone_comp
allcomp-=set(controls.values())
print(f"total compensator candidates (1-hop): {len(allcomp)}")
# how many are currently zero (fresh slack) vs determined
zero=sum(1 for v in allcomp if val[v]==0)
print(f"  of which currently zero (fresh slack): {zero}, nonzero: {len(allcomp)-zero}")
json.dump({'controls':list(controls.values()),'xcone':xcone,'compensators':sorted(allcomp)},
          open('agentC_probe.json','w'))
print("saved agentC_probe.json")
