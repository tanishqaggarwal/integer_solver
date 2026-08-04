#!/usr/bin/env python3
"""Experiment: nudge two handles to force S,T ==0 mod p; measure wiring damage; then try pinning
those handles + re-forward to see if wiring re-heals while residues stay 0."""
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
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best: val[v]=best[v]
forward(); ns['v']=val
F0=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"baseline: {len(lines)-len(F0)}/{len(lines)} ({len(F0)} fail)")
def inv(a): return pow(a%p,p-2,p)
def jac_of(w):
    old=val[w]; val[w]=old+1; forward()
    dS=(val[35389]-Sc)%p; dT=(val[6671]-Tc)%p; val[w]=old; forward(); return dS,dT
Sc=val[35389]; Tc=val[6671]
# pick handles: prefer non-pinned, non-override, with independent (dS,dT)
cone=sorted((anc[35389]|anc[6671])&freeinp)
cand_h=[w for w in cone if not pinned[w]]
# gather jac
J={}
for w in cand_h:
    dS,dT=jac_of(w)
    if dS or dT: J[w]=(dS,dT)
print(f"movable non-pinned handles: {len(J)}")
# find independent pair
hs=list(J); h1=None;h2=None
for i in range(len(hs)):
    for j in range(i+1,len(hs)):
        a,b=J[hs[i]]; c,d=J[hs[j]]
        if (a*d-b*c)%p!=0: h1,h2=hs[i],hs[j]; break
    if h1: break
print(f"chosen handles: h1=x_{h1} {tuple(hex(z)[:10] for z in J[h1])}, h2=x_{h2} {tuple(hex(z)[:10] for z in J[h2])}")
# solve [dS1 dS2; dT1 dT2] [d1;d2] = -(Sc,Tc) mod p
a,b=J[h1]; c,d=J[h2]
det=(a*c- 0)  # placeholder
# matrix columns are handle vectors: [[a, c],[b, d]] * [d1,d2]^T = [-Sc,-Tc]
D=(a*d-b*c)%p; Di=inv(D)
rS=(-Sc)%p; rT=(-Tc)%p
d1=((d*rS - c*rT)*Di)%p
d2=((a*rT - b*rS)*Di)%p
val[h1]=(val[h1]+d1); val[h2]=(val[h2]+d2)
forward(); ns['v']=val
print(f"after nudge: S mod p={val[35389]%p}, T mod p={val[6671]%p}")
F1=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"after nudge (no re-forward-construct): {len(lines)-len(F1)}/{len(lines)} ({len(F1)} fail)")
# check loads mod p
print(f"L1 mod p={val[11150]%p}, L2 mod p={val[25739]%p}, L3 mod p={val[37758]%p}")
