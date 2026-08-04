#!/usr/bin/env python3
"""wire=1 from 39013 with CORRECT wire=1 quotients: x_30317=-L1, x_2936=537773*L3, x_5146=L2/6672769.
Keeps x_16742,x_12186 at baseline -> 12 loads intact. Reports remaining fails (should be 13 unpackings
+ M2-related if L2 not divisible). Then this is the clean base for the final heal."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977
A=load_atoms()
par=list(range(NVARS)); sgn=[1]*NVARS
def find2(x):
    s=1; r=x
    while par[r]!=r: s*=sgn[r]; r=par[r]
    return r,s
def union(a,b,rel):
    ra,sa=find2(a); rb,sb=find2(b)
    if ra==rb: return
    par[ra]=rb; sgn[ra]=rel*sb*sa
for pp in A:
    vs=atom_vars(pp)
    if len(vs)==2 and pp.get((),0)==0:
        v1,v2=sorted(vs); c1=pp.get((v1,),0); c2=pp.get((v2,),0)
        qok=all(pp.get(k,0)==0 for k in pp if isinstance(k,tuple) and len(k)==2)
        if qok and c1!=0 and c2!=0 and abs(c1)==abs(c2):
            rel=(-c2)//c1
            if rel in (1,-1): union(v1,v2,rel)
r0,_=find2(26064); wire={v:find2(v)[1] for v in range(NVARS) if find2(v)[0]==r0}
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
gate_out=set(t for t,_,_ in gates); freeinp=set(v for v in range(NVARS) if v not in gate_out)
best={int(k[2:]):v for k,v in json.load(open('best/new_instance_partial_39013.json')).items()}
val=[0]*NVARS
for k,v in best.items():
    if k<NVARS: val[k]=v
for v,s in wire.items(): val[v]=s*1
# topo order (wire excluded from targets since pinned)
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or v in freeinp or v in wire: ready[v]=True
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
    if ready[t] or t in wire: continue
    definer[t]=gi; order.append(t); ready[t]=True
    for gj in using[t]:
        gu[gj]-=1
        if gu[gj]==0: q.append(gj)
VAR=re.compile(r'x_(\d+)')
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
eqcode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
ns={'__builtins__':{},'v':val}
def forward():
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
forward()
L1=val[11150]; L2=val[25739]; L3=val[37758]
# CORRECT wire=1 quotients
val[30317]=-L1
val[2936]=537773*L3
if L2%6672769==0: val[5146]=L2//6672769; m2ok=True
else: val[5146]=L2//6672769; m2ok=False
forward()
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"wire=1 + correct quotients: {len(lines)-len(F)}/{len(lines)} ({len(F)} fail)")
print(f"  L2%6672769={L2%6672769} (M2 {'OK' if m2ok else 'FAIL'})")
print(f"  core fails: {sorted(i for i in F if i in CORE)}")
print(f"  noncore fails: {sorted(i for i in F if i not in CORE)}")
# check the 12 loads still satisfied?
loads12=[4833,4944,5348,9344,10406,11574,12321,19708,20927,21972,27514,38014]
print(f"  12 loads still OK: {all(i not in F for i in loads12)}")
if True:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('wire1correct.json','w')); print("saved wire1correct.json")
