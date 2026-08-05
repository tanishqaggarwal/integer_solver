#!/usr/bin/env python3
"""For the 26 failing equations, map variable coupling to satisfied equations."""
import json, re
from collections import defaultdict
from propagate import NVARS
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
cand={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39007.json')).items()}
v=[0]*NVARS
for k,x in cand.items():
    if k<NVARS: v[k]=x
codes=[compile(re.sub(r'x_(\d+)', r'v[\1]', L.rsplit('=',1)[0]), '<e>','eval') for L in lines]
ns={'v':v,'__builtins__':{}}
fail=set(); 
eqvars=[]
VAR=re.compile(r'x_(\d+)')
for i,L in enumerate(lines):
    ids=set(int(m) for m in VAR.findall(L))
    eqvars.append(ids)
    if eval(codes[i],ns)!=0: fail.add(i)
print(f"{len(fail)} failing eqs")
# var -> list of eqs containing it
var2eq=defaultdict(list)
for i,ids in enumerate(eqvars):
    for z in ids: var2eq[z].append(i)
# union of vars in failing eqs
failvars=set()
for i in fail: failvars|=eqvars[i]
print(f"failing eqs touch {len(failvars)} distinct vars")
# classify each failing var: appears only in failing eqs (SAFE) vs also in satisfied eqs
safe=[]; coupled=[]
for z in failvars:
    others=[e for e in var2eq[z] if e not in fail]
    if not others: safe.append(z)
    else: coupled.append((z, len(others)))
print(f"SAFE vars (only in failing eqs): {len(safe)}")
print(f"COUPLED vars (also in satisfied eqs): {len(coupled)}")
# how many coupled vars are only lightly coupled (in 1-2 satisfied eqs)?
light=[z for z,n in coupled if n<=2]
print(f"  lightly coupled (<=2 satisfied eqs): {len(light)}")
# per failing equation: how many SAFE vars does it have?
safeset=set(safe)
print("\nPer failing equation: (#vars, #safe-vars, current residual bit-length)")
for i in sorted(fail):
    r=eval(codes[i],ns)
    ns_safe=len(eqvars[i]&safeset)
    print(f"  eq[{i}]: nvars={len(eqvars[i])} safe={ns_safe} resid~2^{r.bit_length()-1 if r else 0}")
json.dump({'safe':sorted(safe),'fail':sorted(fail)}, open('coupling.json','w'))
