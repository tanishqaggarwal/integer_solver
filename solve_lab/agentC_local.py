#!/usr/bin/env python3
"""Examine the equations x_12186 and x_16742 directly touch (their downstream footprint in
satisfied eqs), the compensators available, and each compensator's GLOBAL multi-role degree
(how many eqs it appears in). Identify fresh slack (degree small) for a clean local ℤ solve."""
import json
from agentC_common import (p, gates, order, definer, forward, val, freeinp, ns, lines, eqcode,
                           eqvars, load_best, CORE, downstream_ks, partial_forward, rootcode_of)
from collections import defaultdict

best=load_best(); forward()
gate_defs={t:(rhs,vids) for t,rhs,vids in gates}
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
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
F0=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
deg=defaultdict(int)
for i in range(len(lines)):
    for v in eqvars[i]&freeinp: deg[v]+=1

def footprint(h):
    aff=set(eqbyvar.get(h,()))
    for k in downstream_ks(h): aff|=eqbyvar.get(order[k],set())
    o=val[h]; val[h]=o+1; partial_forward(downstream_ks(h))
    broke=[i for i in aff if i not in F0 and eval(rootcode_of(i),ns)%p!=0]
    val[h]=o; partial_forward(downstream_ks(h))
    return broke

for ctrl in (12186,16742):
    broke=footprint(ctrl)
    print(f"\n===== x_{ctrl}: breaks {len(broke)} sat-eqs when perturbed =====")
    print(f"broken eqs: {sorted(broke)}")
    # for each broken eq, list its free inputs (excluding deep cone) with global degree, and whether
    # eq contains x_ctrl directly or via gate
    deepcone=set()
    for d in [3558,29322,33469,27713,1326,35389,6671]: deepcone|=freecone(d)
    for i in sorted(broke)[:20]:
        fis=eqvars[i]&freeinp
        comp=[(v,deg[v]) for v in fis if v not in deepcone]
        comp.sort(key=lambda x:x[1])
        fresh=[v for v,dd in comp if dd<=2]
        print(f"  eq {i}: |freeinp|={len(fis)}, direct={ctrl in eqvars[i]}, "
              f"fresh(deg<=2) comps={fresh[:6]}, low-deg comps={[(v,dd) for v,dd in comp[:5]]}")
