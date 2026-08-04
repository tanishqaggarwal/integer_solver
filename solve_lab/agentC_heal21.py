#!/usr/bin/env python3
"""Heal the final 12 failures in best_agentA_39021.json (raw assignment). Need x_33462=CONST1,
x_22152=CONST2. Diagnose: are they free/gate? what eqs contain them? what breaks if set?"""
import json, re, sys
from collections import defaultdict
p=2**256-2**32-977
NVARS=38748
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
VAR=re.compile(r'x_(\d+)')
lines=[L.strip() for L in open('../EQUATIONS.txt') if L.strip()]
codes=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
varsets=[frozenset(int(m) for m in VAR.findall(L.rsplit('=',1)[0])) for L in lines]
eqbyvar=defaultdict(set)
for i,vs in enumerate(varsets):
    for v in vs: eqbyvar[v].add(i)
# gate outputs
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for L in f:
        gate_out.add(json.loads(L)['t'])
a21={int(k[2:]):v for k,v in json.load(open('best_agentA_39021.json')).items()}
V=[0]*NVARS
for k,v in a21.items(): V[k]=v
ns={'v':V,'__builtins__':{}}
def failing():
    return [i for i in range(len(lines)) if eval(codes[i],ns)!=0]
F0=failing()
print(f"agentA_39021: {len(lines)-len(F0)}/{len(lines)} ({len(F0)} fail): {sorted(F0)}")
print(f"x_33462: free={33462 not in gate_out}, current={V[33462]}, target CONST1")
print(f"  matches CONST1 already: {V[33462]==CONST1}; #eqs={len(eqbyvar[33462])}")
print(f"x_22152: free={22152 not in gate_out}, current={V[22152]}, target CONST2")
print(f"  matches CONST2 already: {V[22152]==CONST2}; #eqs={len(eqbyvar[22152])}")
# do the 12 failing contain x_33462/x_22152?
for i in F0:
    print(f"  eq {i}: has_33462={33462 in varsets[i]} has_22152={22152 in varsets[i]} nvars={len(varsets[i])}")
# set the two loads, see what breaks
V[33462]=CONST1; V[22152]=CONST2
F1=failing()
print(f"\nafter setting x_33462=CONST1,x_22152=CONST2: {len(lines)-len(F1)}/{len(lines)} ({len(F1)} fail)")
fixed=set(F0)-set(F1); broke=set(F1)-set(F0)
print(f"  fixed: {sorted(fixed)}")
print(f"  newly broken ({len(broke)}): {sorted(broke)}")
