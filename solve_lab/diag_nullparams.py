#!/usr/bin/env python3
"""Run forward_construct's logic; find UNDETERMINED free inputs (null params). Test which move
x_3558/x_29322 residues AND how many equations they break (0 broken => free null direction)."""
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
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
eqvars=[set(int(m) for m in VAR.findall(L)) for L in lines]
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best: val[v]=best[v]
forward(); ns['v']=val
# determined = free inputs with nonzero best value OR pinned; undetermined = value 0 and not pinned
undet=[v for v in freeinp if not pinned[v] and val[v]==0]
det=[v for v in freeinp if not pinned[v] and val[v]!=0]
print(f"free inputs: {len(freeinp)}; pinned {sum(1 for v in freeinp if pinned[v])}; determined(nonzero) {len(det)}; undetermined(zero) {len(undet)}")
F0=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
# test each undetermined: does it move x_3558 or x_29322 residue? does it break eqs?
movers=[]
for v in undet:
    val[v]=1; forward(); ns['v']=val
    d3=(val[3558]-0)%p if val[3558]%p else 0  # residue change proxy: recompute
    # measure residue change: compare x_3558 mod p vs baseline (baseline x_3558%p known)
    val[v]=0
forward(); ns['v']=val
b3=val[3558]%p; b29=val[29322]%p
for v in undet:
    val[v]=1; forward(); ns['v']=val
    d3=(val[3558]%p-b3)%p; d29=(val[29322]%p-b29)%p
    if d3 or d29:
        Fb=sum(1 for i in range(len(lines)) if i not in F0 and eval(eqcode[i],ns)!=0)
        movers.append((v,d3!=0,d29!=0,Fb))
    val[v]=0
forward()
print(f"undetermined movers of x_3558/x_29322 residue: {len(movers)}")
for v,m3,m29,fb in movers[:20]:
    print(f"  x_{v}: moves3558={m3}, moves29322={m29}, breaks={fb} eqs")
