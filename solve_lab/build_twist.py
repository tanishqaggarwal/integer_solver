#!/usr/bin/env python3
"""Twist solution via MUX activation:
control bits (x_7715,x_34554)=(1,0) select quadrant x_34606; MUX routes free inputs directly:
x_37892=x_16742=C2, x_13682=x_12186=C1. Activate x_7715 by one free-input flip; keep x_34554=0.
Also pin square-only vars x_24468=C1, x_18956=C2."""
import json, re, ast, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS

hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
ACT=int(sys.argv[1]) if len(sys.argv)>1 else 47
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
# pins from single-var atoms
val=[0]*NVARS; pinned=[False]*NVARS
for p in A:
    vs=atom_vars(p)
    if len(vs)==1:
        v=next(iter(vs)); c0=p.get((),0); c1=p.get((v,),0); c2=p.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]:
            val[v]=(-c0)//c1; pinned[v]=True
# OVERRIDES: activation + MUX routing + square-only huge pins
override={ACT:1, 16742:C2, 12186:C1, 24468:C1, 18956:C2}
for v,x in override.items(): val[v]=x; pinned[v]=True
# greedy topo orientation
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand)
ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or pinned[v]: ready[v]=True
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
code={t:compile(VAR.sub(r'v[\1]',gates[definer[t]][1]),'<r>','eval') for t in order}
ns={'__builtins__':{}}
for t in order:
    ns['v']=val; val[t]=eval(code[t],ns)
out={f"x_{i}":val[i] for i in range(NVARS)}
json.dump(out, open('twist_solved.json','w'))
# verify: control bits + routed outputs
print(f"activator x_{ACT}=1 -> x_7715={val[7715]} x_34554={val[34554]} x_9274={val[9274]} x_15298={val[15298]} x_34606={val[34606]}")
print(f"x_37892={val[37892]==C2} (want C2) x_13682={val[13682]==C1} (want C1)")
# equation-level count
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
codes=[compile(re.sub(r'x_(\d+)', r'v[\1]', L.rsplit('=',1)[0]), '<e>','eval') for L in lines]
ns2={'v':val,'__builtins__':{}}
fails=[i for i,c in enumerate(codes) if eval(c,ns2)!=0]
print(f"EQUATION LEVEL: {len(lines)-len(fails)}/{len(lines)} ({len(fails)} fail): {fails[:40]}")
json.dump(fails, open('twist_fails.json','w'))
