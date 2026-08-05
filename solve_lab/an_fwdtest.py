#!/usr/bin/env python3
"""Forward-reconstruct from agentA's free inputs, but x_4432,x_7068 reverted to 39013 values.
Tests: does core survive + 11 get fixed?"""
import json,re,ast
from collections import defaultdict,deque
p=2**256-2**32-977; NVARS=38748
VAR=re.compile(r'x_(\d+)')
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'],d['rhs'],tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
ready=[False]*NVARS
for v in range(NVARS):
    if v in freeinp: ready[v]=True
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
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[t]][1]),'<r>','eval') for t in order]
def loadd(path):
    d=json.load(open(path)); out={}
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); out[idx]=int(val)
    return out
vA=loadd('best_agentA_39022.json'); v013=loadd('best/new_instance_partial_39013.json')
# count free inputs NOT defined by any gate that agentA covers
print(f"free inputs total={len(freeinp)}, gate outputs={len(order)}, undefined gate-targets={len(gate_out)-len(order)}")
val=[0]*NVARS
for v in freeinp: val[v]=vA.get(v,0)
# revert the two
val[4432]=v013[4432]; val[7068]=v013[7068]
ns={'v':val,'__builtins__':{}}
for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
# check all eqs
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
F=[i for i,c in enumerate(eqcode) if eval(c,ns)!=0]
print(f"forward-reconstruct (x_4432,x_7068 reverted): {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)")
print("fails:",F[:40])
# Now the SAME but WITHOUT reverting (pure agentA frees -> forward): baseline
for v in freeinp: val[v]=vA.get(v,0)
ns={'v':val,'__builtins__':{}}
for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
F2=[i for i,c in enumerate(eqcode) if eval(c,ns)!=0]
print(f"\nbaseline (agentA frees, pure forward): {len(lines)-len(F2)}/{len(lines)} ({len(F2)} fail)")
print("fails:",F2[:40])
