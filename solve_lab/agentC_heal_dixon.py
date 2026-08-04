#!/usr/bin/env python3
"""Single-shot rank-Dixon heal of the final 12 in best_agentA_39021, using agentA_harness (whose
forward() matches agentA's construction). Set loads to CONST, forward; heal the resulting breaks
with ONE Dixon solve over handles that don't feed the load cone (no feedback -> loads stay fixed)."""
import json, re, sys
from collections import defaultdict, deque
from agentA_harness import (p, order, definer, gates, freeinp, backward_cone, load_solution,
                            forward, eqcode, eqvars, lines, NEQ, NVARS)
sys.setrecursionlimit(1000000)
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
CORE={2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892}
QUOT={30317,2936,5146}
gvids={t:gates[definer[t]][2] for t in order}
ns={'__builtins__':{}}
def count(v):
    ns['v']=v; return set(i for i in range(NEQ) if eval(eqcode[i],ns)!=0)
def set_quot(v):
    if v[11150]%p==0: v[30317]=-(v[11150])//p
    if (537773*v[37758])%p==0: v[2936]=(537773*v[37758])//p
    if v[25739]%(6672769*p)==0: v[5146]=v[25739]//(6672769*p)
def free_cone(r):
    seen=set(); st=[r]
    while st:
        u=st.pop()
        if u in seen: continue
        seen.add(u)
        for w in gvids.get(u,()):
            if w not in seen: st.append(w)
    return set(u for u in seen if u in freeinp)

base=load_solution('best/new_instance_partial_39013.json'); forward(base)
v=load_solution('best_agentA_39021.json')
before=count(v)
print(f"agentA_39021 as-loaded: {NEQ-len(before)}/{NEQ} ({len(before)} fail)")
v[22152]=CONST2; v[33462]=CONST1; forward(v); set_quot(v)
cur=count(v)
print(f"after loads=CONST + forward + set_quot: {NEQ-len(cur)}/{NEQ} ({len(cur)} fail) core={len(cur&CORE)}")
print(f"broken: {sorted(cur)}")
# load cone: free inputs feeding the (now fixed) 12 load eqs + core + the two loads
LOAD12=[4833,4944,5348,9344,10406,11574,12321,19708,20927,21972,27514,38014]
loadcone=set()
for i in LOAD12:
    for var in eqvars[i]:
        loadcone|=free_cone(var);
        if var in freeinp: loadcone.add(var)
corecone=free_cone(35389)|free_cone(6671)|{35389,6671}
protect=loadcone|corecone|QUOT|{22152,33462,12186,14853,16742,24908}
print(f"load-cone frees: {len(loadcone)}, core-cone frees: {len(corecone)}, protect total: {len(protect)}")
# candidate handles: free inputs feeding the broken eqs, minus protect
Hset=set()
for i in cur:
    for var in eqvars[i]:
        Hset|=free_cone(var)
        if var in freeinp: Hset.add(var)
H=sorted(Hset-protect)
print(f"candidate non-load-cone handles: {len(H)}")
json.dump({'broken':sorted(cur),'handles':H,'loadcone':sorted(loadcone)}, open('agentC_healdix.json','w'))
