#!/usr/bin/env python3
"""Reusable forward-reconstruct + tangent-linear Jacobian harness on the wire=p (agentA) branch.
Free inputs -> forward() -> all gate outputs. Core auto-maintained. Heal target eqs by moving
free inputs, GF(p) least-effort solve, iterate. Pin specified free inputs."""
import json,re,sys
from collections import defaultdict,deque
p=2**256-2**32-977; NVARS=38748
VAR=re.compile(r'x_(\d+)')
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'],d['rhs'],tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
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
anc=defaultdict(set)
for v in freeinp: anc[v]={v}
for t in order:
    _,rhs,vids=gates[definer[t]]
    s=set()
    for u in vids: s|=anc[u]
    anc[t]=s
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[frozenset(int(m) for m in VAR.findall(L)) for L in lines]
def loadd(path):
    d=json.load(open(path)); out={}
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); out[idx]=int(val)
    return out
val=[0]*NVARS
def forward():
    ns={'v':val,'__builtins__':{}}
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
def fails():
    ns={'v':val,'__builtins__':{}}
    return [i for i,c in enumerate(eqcode) if eval(c,ns)!=0]
def resid(idxs):
    ns={'v':val,'__builtins__':{}}
    return {i:eval(eqcode[i],ns) for i in idxs}

if __name__=='__main__':
    vA=loadd('best_agentA_39022.json'); v013=loadd('best/new_instance_partial_39013.json')
    for v in freeinp: val[v]=vA.get(v,0)
    forward()
    F=fails()
    print(f"agentA baseline: {len(F)} fail: {F}")
    # descendants of x_4432,x_7068
    desc=set()
    for t in order:
        if (anc[t]&{4432,7068}): desc.add(t)
    print(f"gate outputs downstream of x_4432/x_7068: {len(desc)}")
