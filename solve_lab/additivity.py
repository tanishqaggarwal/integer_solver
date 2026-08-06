#!/usr/bin/env python3
"""Test whether bit loads are additive (linear message) or bilinear.
For monitor eqs (those touching any 256 bit), compare resid({A,B}) vs resid({A})+resid({B})-resid({})."""
import json, re
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
relbits=json.load(open('boolbits.json'))['relbits']  # 258 spine bits
# use x_7715 cone bits as the activators
bits7=json.load(open('act7715.json'))['free7']
gate_out=set(t for t,_,_ in gates)
for v in range(NVARS):
    if v not in gate_out: pinned[v]=True
# topo
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
gcode=[compile(VAR.sub(r'v[\1]',gates[definer[order[k]]][1]),'<r>','eval') for k in range(len(order))]
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
# monitor eqs: those touching any relbit
bitset=set(relbits)|set(bits7)
mon=[i for i,L in enumerate(lines) if set(int(m) for m in VAR.findall(L)) & bitset]
print(f"monitor eqs: {len(mon)}")
mcode=[compile(VAR.sub(r'v[\1]',lines[i].rsplit('=',1)[0]),'<e>','eval') for i in mon]
ns={'__builtins__':{}}
def resid(setbits):
    val=val0[:]
    for b in setbits: val[b]=1
    ns['v']=val
    for k,t in enumerate(order): val[t]=eval(gcode[k],ns)
    return [eval(c,ns) for c in mcode]
r0=resid([])
# pick two activators
Abit, Bbit = bits7[0], bits7[1]
rA=resid([Abit]); rB=resid([Bbit]); rAB=resid([Abit,Bbit])
mism=[]
for k in range(len(mon)):
    pred=rA[k]+rB[k]-r0[k]
    if rAB[k]!=pred: mism.append(mon[k])
print(f"pair (x_{Abit},x_{Bbit}): additivity violations {len(mism)}/{len(mon)}: {mism[:10]}")
# try a few more pairs
import itertools
for Abit,Bbit in [(bits7[2],bits7[3]),(bits7[0],bits7[5]),(bits7[10],bits7[20])]:
    rA=resid([Abit]); rB=resid([Bbit]); rAB=resid([Abit,Bbit])
    mism=[mon[k] for k in range(len(mon)) if rAB[k]!=rA[k]+rB[k]-r0[k]]
    print(f"pair (x_{Abit},x_{Bbit}): violations {len(mism)}: {mism[:8]}")
