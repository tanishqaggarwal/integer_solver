#!/usr/bin/env python3
"""Newton-iterate on handles to drive S,T -> 0 mod p (holding other free inputs). Watch convergence."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
val=[0]*NVARS; pinned=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val[v]=(-c0)//c1; pinned[v]=True
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
override={24601:1, 2081:1, 30213:C2, 22162:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x; pinned[v]=True
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or v in freeinp or pinned[v]: ready[v]=True
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
VAR=re.compile(r'x_(\d+)')
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
anc=defaultdict(set)
for v in freeinp: anc[v]={v}
for k,t in enumerate(order):
    s=set()
    for u in gates[definer[t]][2]: s|=anc[u]
    anc[t]=s
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best: val[v]=best[v]
forward()
def inv(a): return pow(a%p,p-2,p)
# Newton on two handles h1,h2 (deep, independent). Recompute jac each iter.
h1,h2=91,12186
for it in range(25):
    Sc=val[35389]%p; Tc=val[6671]%p
    if Sc==0 and Tc==0: print(f"iter {it}: CONVERGED S=T=0 mod p"); break
    # jac
    o1=val[h1]; val[h1]=o1+1; forward(); a=(val[35389]-Sc)%p; b=(val[6671]-Tc)%p; val[h1]=o1
    o2=val[h2]; val[h2]=o2+1; forward(); c=(val[35389]-Sc)%p; d=(val[6671]-Tc)%p; val[h2]=o2
    forward()
    D=(a*d-b*c)%p
    if D==0: print(f"iter {it}: singular jac"); break
    Di=inv(D); rS=(-Sc)%p; rT=(-Tc)%p
    d1=((d*rS-c*rT)*Di)%p; d2=((a*rT-b*rS)*Di)%p
    val[h1]+=d1; val[h2]+=d2; forward()
    if it<6 or it%5==0: print(f"iter {it}: S mod p bits={ (val[35389]%p).bit_length() }, T mod p bits={ (val[6671]%p).bit_length() }")
print(f"FINAL: S mod p = {val[35389]%p}")
print(f"       T mod p = {val[6671]%p}")
