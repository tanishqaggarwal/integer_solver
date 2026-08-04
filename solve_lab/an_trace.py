#!/usr/bin/env python3
import json,re
from collections import defaultdict,deque
p=2**256-2**32-977; NVARS=38748
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'],d['rhs'],tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
# topo order via forward_construct logic (simplified: pick first definer reachable)
# Build definer by Kahn using ANY gate per target
ready=[False]*NVARS
# free inputs + pinned single-var are ready; approximate: freeinp ready
override={24601,2081,30213,22162,24468,18956}
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
anc=defaultdict(set)
for v in freeinp: anc[v]={v}
for t in order:
    _,rhs,vids=gates[definer[t]]
    s=set()
    for u in vids: s|=anc[u]
    anc[t]=s
def loadd(path):
    d=json.load(open(path)); out={}
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); out[idx]=int(val)
    return out
vA=loadd('best_agentA_39022.json'); v013=loadd('best/new_instance_partial_39013.json')
changed=set(k for k in set(vA)|set(v013) if vA.get(k,0)!=v013.get(k,0))
changed_free=changed&freeinp
print(f"total changed={len(changed)}, changed free inputs={len(changed_free)}: {sorted(changed_free)}")
for tgt in [4432,7068]:
    isfree = tgt in freeinp
    print(f"\nx_{tgt}: freeinput={isfree}, #free-ancestors={len(anc.get(tgt,set()))}")
    anc_t=anc.get(tgt,set())
    conflict=anc_t & changed_free
    print(f"  free ancestors changed by agentA: {sorted(conflict)}")
