#!/usr/bin/env python3
"""Heal agentA_39021: set x_33462=CONST1, x_22152=CONST2, then targeted forward-propagate ONLY
their downstream gate cone (keep all other agentA values). Check + save if >39021."""
import json, re, sys
import agentC_common as AC
from agentC_common import (p, order, gcode, definer, gates, downstream_ks, val, ns, lines, eqcode,
                           NVARS, CORE, freeinp, posof)
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
a21={int(k[2:]):v for k,v in json.load(open('best_agentA_39021.json')).items()}
for i in range(NVARS): val[i]=a21.get(i,0)
ns['v']=val
F0=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"start: {len(lines)-len(F0)}/{len(lines)} fail={sorted(F0)}")
# is core currently satisfied? (S,T)
print(f"S=x_35389={val[35389]%p}, T=x_6671={val[6671]%p}")
# set the two loads
val[33462]=CONST1; val[22152]=CONST2
# downstream gate cone (topo positions) of the two changed free inputs
D=sorted(set(downstream_ks(33462)) | set(downstream_ks(22152)))
print(f"downstream gates to repropagate: {len(D)}")
# check whether core gates (35389,6671) are in the cone
core_gates_in={t for t in (35389,6671,33469,29322,3558,27713,1326) if t in posof and posof[t] in set(D)}
print(f"core gates in propagation cone (should be empty): {core_gates_in}")
for k in D:
    val[order[k]]=eval(gcode[k], ns)
ns['v']=val
F1=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
core=[i for i in F1 if i in CORE]; nc=[i for i in F1 if i not in CORE]
print(f"after set+propagate: {len(lines)-len(F1)}/{len(lines)} ({len(F1)} fail)")
print(f"  core-fail={len(core)}: {sorted(core)}")
print(f"  noncore-fail={len(nc)}: {sorted(nc)}")
print(f"S=x_35389={val[35389]%p}, T=x_6671={val[6671]%p}")
if len(F1)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS) if val[i]!=0}, open('best_agentC_39033.json','w'))
    print("*** FULL WIN — saved best_agentC_39033.json ***")
elif len(F1)<12:
    json.dump({f"x_{i}":val[i] for i in range(NVARS) if val[i]!=0}, open(f'best_agentC_{len(lines)-len(F1)}.json','w'))
    print(f"improved -> saved best_agentC_{len(lines)-len(F1)}.json")
