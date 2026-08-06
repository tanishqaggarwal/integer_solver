#!/usr/bin/env python3
"""Wire=1 base (39013+wire=1, core via quotients, loads intact, S,T NOT zeroed so x_14853 untouched
-> NO 11-fail wall). Fix M2 by shifting x_3558 to root mod 6672769 via x_31339 (x_24908 knob), NOT
x_16742/x_14853. Re-set quotients. Report remaining fails (should be 13 unpackings + x_31339 ripple)."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
sys.setrecursionlimit(1000000)
p=2**256-2**32-977; M=6672769
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
def setquot():
    val[30317]=-val[11150]; val[2936]=537773*val[37758]
    if val[25739]%M==0: val[5146]=val[25739]//M
forward(); setquot(); forward()
CORE=set([2071,4573,7123,7469,11854,13660,15299,16622,17726,21382,22093,25480,25539,28653,29437,31061,32894,32916,34517,34892])
F=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"before M2 shift: {len(F)} fail (core {len([i for i in F if i in CORE])}); L2%M={val[25739]%M}, y=x_3558%M={val[3558]%M}")
# knob coeff of x_31339 in x_3558
o=val[31339]; val[31339]=o+1; forward(); coeff=(val[3558]-val[3558]) # placeholder
val[31339]=o; forward()
b3558=val[3558]
val[31339]=o+1; forward(); coeff=(val[3558]-b3558); val[31339]=o; forward()
print(f"x_31339 -> x_3558 coeff = {coeff}")
# solve L2%M=0: 6506 + 777865*y + 3186439*y^2 = 0 mod M. roots y in {2783706,5883594}
cur_y=val[3558]%M
for root in [2783706,5883594]:
    dY=(root-cur_y)%M
    # need x_3558 += dY  => x_31339 += dY * inv(coeff) mod M
    if coeff%M==0: continue
    dk=(dY*pow(coeff%M,M-2,M))%M
    val[31339]=o+dk; forward(); setquot(); forward()
    l2m=val[25739]%M
    F2=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    print(f"root {root}: dk={dk}, L2%M={l2m}, total fail={len(F2)} (core {len([i for i in F2 if i in CORE])})")
    if l2m==0:
        nc=[i for i in F2 if i not in CORE]
        print(f"  M2 FIXED! remaining: core {len([i for i in F2 if i in CORE])}, noncore {len(nc)}: {sorted(nc)}")
        json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('wire1_m2fixed.json','w')); print("  saved wire1_m2fixed.json")
        break
    val[31339]=o; forward(); setquot(); forward()
