#!/usr/bin/env python3
"""The 12 fails in 39021 reduce to x_24601*(x_33462-CONST1)=0 and x_24601*(x_22152-CONST2)=0.
Extract CONST1,CONST2; set x_33462,x_22152 to zero these; forward-eval; verify. Maybe FULL SOLVE."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
s={int(k[2:]):v for k,v in json.load(open('best_agentA_39021.json')).items()}
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
VAR=re.compile(r'x_(\d+)')
# extract CONST1 (with x_33462) and CONST2 (with x_22152) from any gadget
def get_const(pv):
    for i in [4833,4944,5348,9344,10406,11574,12321,19708,20927,21972,27514,38014]:
        for m in re.finditer(r'x_%d\)\s*-\s*\(?(\d+)'%pv, lines[i]):
            return int(m.group(1))
        # try (x_pv - CONST) pattern
    for i in [4833,4944,5348,9344,10406,11574,12321,19708,20927,21972,27514,38014]:
        mm=re.search(r'\(x_%d\)-\((\d+)\)'%pv, lines[i].replace(' ',''))
        if mm: return int(mm.group(1))
    return None
C1=get_const(33462); C2=get_const(22152)
print(f"x_24601={s.get(24601)}, x_33462={s.get(33462)} (bits {s.get(33462,0).bit_length()}), x_22152={s.get(22152)}")
print(f"CONST1={C1} (bits {C1.bit_length() if C1 else 0})")
print(f"CONST2={C2} (bits {C2.bit_length() if C2 else 0})")
# rebuild full forward-eval to be safe (x_33462,x_22152 are free inputs feeding gates)
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
val=[0]*NVARS
for k,v in s.items():
    if k<NVARS: val[k]=v
# topo order of gates
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or v in freeinp: ready[v]=True
gu=[0]*len(gates); using=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates):
    u=0
    for v in vids:
        if not ready[v]: u+=1
        using[v].append(gi)
    gu[gi]=u
definer={}; order=[]
q=deque(gi for gi in range(len(gates)) if gu[gi]==0)
while q:
    gi=q.popleft(); t,rhs,vids=gates[gi]
    if ready[t]: continue
    definer[t]=gi; order.append(t); ready[t]=True
    for gj in using[t]:
        gu[gj]-=1
        if gu[gj]==0: q.append(gj)
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
ns={'__builtins__':{},'v':val}
def forward():
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
forward()
F0=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"before: {len(lines)-len(F0)}/{len(lines)} ({len(F0)} fail)")
# set free inputs
if C1 is not None: val[33462]=C1
if C2 is not None: val[22152]=C2
forward()
F1=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"after set x_33462=CONST1, x_22152=CONST2: {len(lines)-len(F1)}/{len(lines)} ({len(F1)} fail): {sorted(F1)[:20]}")
if len(F1)==0:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('FULL_SOLVED.json','w')); print("*** FULL SOLUTION — saved FULL_SOLVED.json ***")
elif len(F1)<12:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('fix12_partial.json','w')); print(f"improved to {len(lines)-len(F1)}")
