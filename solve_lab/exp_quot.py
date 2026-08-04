#!/usr/bin/env python3
"""After minimal residue fix (S=T=0 mod p): check load divisibilities, set quotient handles
x_30317,x_5146,x_2936 to zero M1,M2,M3. Report which core eqs clear + M2 divisibility status."""
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
ns={'__builtins__':{}}
def forward():
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
for v in freeinp:
    if v in best: val[v]=best[v]
forward()
val[14853]=val[14853]-(val[29322]%p); forward()
val[16742]=val[16742]+(val[3558]%p); forward()
L1=val[11150]; L2=val[25739]; L3=val[37758]
print(f"L1%p={L1%p}, L2%p={L2%p}, L3%p={L3%p}")
print(f"L1/p integer: {L1%p==0}; L3: 537773*L3 % p == {(537773*L3)%p}")
print(f"L2 % (6672769*p) == {L2 % (6672769*p)}  (need 0 for integer x_5146)")
print(f"L2/p % 6672769 == {(L2//p) % 6672769}")
# set quotient handles
val[30317]=-(L1//p)     # M1: x_11150 + p*x_30317 = 0
val[2936]=(537773*L3)//p  # M3: 537773*L3 - p*x_2936 = 0
if L2 % (6672769*p)==0:
    val[5146]=L2//(6672769*p)  # M2
    print("x_5146 set exactly")
else:
    val[5146]=L2//(6672769*p)  # floor, M2 won't be exactly 0
    print("x_5146 floored (M2 residual)")
forward(); ns['v']=val
F=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
print(f"total: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)")
print(f"core still failing: {sorted(F&CORE)}")
print(f"non-core failing: {sorted(F-CORE)}")
