#!/usr/bin/env python3
"""Rebuild forward-construction. Confirm loads = combos of S=x_35389,T=x_6671. Compute S,T mod p.
Find free inputs in backward cones of S,T; probe Jacobian mod p (which free inputs steer S,T residues)."""
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
# load best partial into free inputs
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best: val[v]=best[v]
forward()
S=val[35389]; T=val[6671]
print(f"S=x_35389 mod p = {S%p}  (val bits {S.bit_length()})")
print(f"T=x_6671  mod p = {T%p}  (val bits {T.bit_length()})")
print(f"L1=x_11150 check: {val[11150]==8646263*S+1073965*T},  L1 mod p={val[11150]%p}")
print(f"L2=x_25739 check: {val[25739]==10159099*S+6926539*T}, L2 mod p={val[25739]%p}")
print(f"L3=x_37758 check: {val[37758]==8272701*S+5921311*T},  L3 mod p={val[37758]%p}")
# minors
import itertools
rows=[(8646263,1073965),(10159099,6926539),(8272701,5921311)]
for (r1,r2) in itertools.combinations(rows,2):
    d=(r1[0]*r2[1]-r1[1]*r2[0])%p
    print(f"  minor det {r1},{r2} mod p = {d} ({'INVERTIBLE' if d!=0 else 'ZERO'})")
# free inputs in cones of S,T
coneS=anc[35389]&freeinp; coneT=anc[6671]&freeinp
print(f"\n|cone(S) free|={len(coneS)}, |cone(T) free|={len(coneT)}, shared={len(coneS&coneT)}")
# Jacobian mod p: perturb each free input in cone by 1, measure dS mod p, dT mod p
S0=S%p; T0=T%p
jac=[]
for w in sorted(coneS|coneT):
    old=val[w]; val[w]=old+1; forward()
    dS=(val[35389]-S)%p; dT=(val[6671]-T)%p
    val[w]=old
    if dS or dT: jac.append((w,dS,dT))
forward()
print(f"free inputs that move (S mod p, T mod p): {len(jac)}")
for w,dS,dT in jac[:30]:
    print(f"   x_{w}: dS={dS if dS<10**6 else hex(dS)[:12]+'..'}, dT={dT if dT<10**6 else hex(dT)[:12]+'..'}")
