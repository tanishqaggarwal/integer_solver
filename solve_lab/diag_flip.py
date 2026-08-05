#!/usr/bin/env python3
"""Quick feasibility read: flip each of the 78 S,T-cone boolean bits (0<->1) from the best solution,
forward-eval, count satisfied + measure S,T residue change. How many flips are 'free' (keep >=39000)?"""
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
gate_defs={t:vids for t,rhs,vids in gates}
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
val=[0]*NVARS; pinned=[False]*NVARS
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==1:
        v=next(iter(vs)); c0=pp.get((),0); c1=pp.get((v,),0); c2=pp.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val[v]=(-c0)//c1; pinned[v]=True
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
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
def fcone(root):
    seen=set(); lv=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x]: st.append(u)
        elif x in freeinp: lv.add(x)
    return lv
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best: val[v]=best[v]
forward(); ns['v']=val
base_ok=sum(1 for i in range(len(lines)) if eval(eqcode[i],ns)==0)
S0=val[35389]%p; T0=val[6671]%p
bb=set(json.load(open('boolbits.json'))['boolvars'])
cone=(fcone(35389)|fcone(6671))&bb&freeinp
cone=sorted(cone)
print(f"baseline satisfied={base_ok}, S%p!=0={S0!=0}, T%p!=0={T0!=0}; testing {len(cone)} cone-boolean bits")
freeflips=[]; dS_list={}
for v in cone:
    old=val[v]; val[v]=1-old if old in (0,1) else (0 if old else 1)
    forward(); ns['v']=val
    ok=sum(1 for i in range(len(lines)) if eval(eqcode[i],ns)==0)
    dS=(val[35389]%p-S0)%p; dT=(val[6671]%p-T0)%p
    val[v]=old
    if ok>=39000: freeflips.append((v,ok,dS!=0,dT!=0))
forward()
print(f"bits whose single flip keeps >=39000 satisfied: {len(freeflips)}")
for v,ok,ds,dt in freeflips[:20]: print(f"  x_{v}: ok={ok}, movesS={ds}, movesT={dt}")
