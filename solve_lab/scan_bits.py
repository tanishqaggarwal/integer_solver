#!/usr/bin/env python3
"""Scan all 178 x_7715-cone activator bits: flip each, forward-eval, count equation failures.
Find minimal-collateral activation."""
import json, re, sys
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
val0=[0]*NVARS; pinned=[False]*NVARS
for p in A:
    vs=atom_vars(p)
    if len(vs)==1:
        v=next(iter(vs)); c0=p.get((),0); c1=p.get((v,),0); c2=p.get((v,v),0)
        if c2==0 and c1!=0 and (-c0)%c1==0 and not pinned[v]: val0[v]=(-c0)//c1; pinned[v]=True
bits=json.load(open('act7715.json'))['free7']
for b in bits: pinned[b]=True
# topo (bits pinned as inputs)
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
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
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[t]][1]),'<r>','eval') for t in order]
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
ecode=[compile(VAR.sub(r'v[\1]',L.rsplit('=',1)[0]),'<e>','eval') for L in lines]
ns={'__builtins__':{}}
def run(setbits):
    val=val0[:]
    for b,x in setbits.items(): val[b]=x
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
    return val
# baseline
val=run({}); ns['v']=val
base_fail=[i for i,c in enumerate(ecode) if eval(c,ns)!=0]
print(f"baseline (all bits 0): {len(lines)-len(base_fail)}/{len(lines)} fail={len(base_fail)}", flush=True)
best=[]
for bi,b in enumerate(bits):
    val=run({b:1}); ns['v']=val
    if val[7715]!=1: continue
    nf=0
    for c in ecode:
        if eval(c,ns)!=0:
            nf+=1
            if nf>len(base_fail)+2: break
    best.append((nf,b))
    if bi%40==0: print(f"  scanned {bi}/178 ...", flush=True)
best.sort()
print(f"\nBest single-bit activations (fail-count, bit):")
for nf,b in best[:15]: print(f"  x_{b}: {nf} fail")
json.dump(best, open('scan_result.json','w'))
