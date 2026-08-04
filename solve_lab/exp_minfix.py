#!/usr/bin/env python3
"""Minimal fix: x_14853 -= (x_29322 mod p) and x_16742 += (x_3558 mod p) so both control gates
become ==0 mod p (changes < p). Measure wiring damage from the best solution; list broken eqs."""
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
F0=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
print(f"baseline: {len(lines)-len(F0)}/{len(lines)}  fail={sorted(F0)}")
# minimal residue fixes
val[14853]=val[14853]-(val[29322]%p)
forward()
val[16742]=val[16742]+(val[3558]%p)
forward(); ns['v']=val
print(f"after fix: x_29322 mod p={val[29322]%p}, x_3558 mod p={val[3558]%p}")
print(f"S mod p={val[35389]%p}, T mod p={val[6671]%p}")
F1=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
print(f"after fix: {len(lines)-len(F1)}/{len(lines)} ({len(F1)} fail)")
newbreak=sorted(F1-F0); fixed=sorted(F0-F1)
print(f"newly broken ({len(newbreak)}): {newbreak}")
print(f"core fixed ({len(fixed)}): {fixed}")
# which free inputs are in the newly-broken eqs (candidates for greedy re-heal)
bh=set()
for i in newbreak: bh|=(eqvars[i]&freeinp)
print(f"free inputs in newly-broken eqs: {len(bh)}")
