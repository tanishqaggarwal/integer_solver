#!/usr/bin/env python3
"""Tie difference-square partners: for each huge square-gate y^2 with y = a - b, set the free-input
side so a == b (difference -> 0), killing the 2^591 monsters. Iterate with forward-eval; verify."""
import json, re, ast
from collections import defaultdict, deque
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
hc=json.load(open('huge_consts.json')); C1=int(hc['C1']); C2=int(hc['C2'])
A=load_atoms()
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gates.append((d['t'], d['rhs'], tuple(d['vids'])))
val=[0]*NVARS
for k,x in json.load(open('best/new_instance_partial_39013.json')).items():
    v=int(k[2:])
    if v<NVARS: val[v]=x
gate_defs=defaultdict(list)
for t,rhs,vids in gates: gate_defs[t].append((rhs,vids))
gate_out=set(gate_defs); freeinp=set(v for v in range(NVARS) if v not in gate_out)
# topo order for forward-eval (single-def)
pinned=[False]*NVARS
for v in freeinp: pinned[v]=True
cand=defaultdict(list)
for gi,(t,rhs,vids) in enumerate(gates): cand[t].append(gi)
targets=set(cand); ready=[False]*NVARS
for v in range(NVARS):
    if v not in targets or v in freeinp: ready[v]=True
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
# free-input ancestors for choosing settable side
anc=defaultdict(set)
for v in freeinp: anc[v]={v}
for k,t in enumerate(order):
    _,vids=gate_defs[t][0] if False else (None, gates[definer[t]][2])
    s=set()
    for u in vids: s|=anc[u]
    anc[t]=s
# find square gates y=a-b
squares={}
for t in gate_out:
    rhs,vids=gate_defs[t][0]
    m=re.match(r'x_(\d+) \* x_(\d+)$', rhs)
    if m and m.group(1)==m.group(2): squares[t]=int(m.group(1))
def diffpair(y):
    rhs,vids=gate_defs.get(y,[('',())])[0]
    m=re.match(r'x_(\d+) - x_(\d+)$', rhs)
    if m: return int(m.group(1)), int(m.group(2))
    return None
forward(); ns['v']=val
Fq=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"start: {len(lines)-len(Fq)}/{len(lines)} ({len(Fq)} fail)", flush=True)
for it in range(40):
    forward(); ns['v']=val
    changed=False
    for t,y in squares.items():
        if abs(val[t])<2**260: continue
        pr=diffpair(y)
        if not pr: continue
        a,b=pr
        if val[a]==val[b]: continue
        # set the free-input side to match the other. Prefer the side whose value is 0 / smaller cone.
        # choose target free input: a single free input in a's or b's cone we can set to zero the diff.
        # simplest: if b is free, set b=val[a]; elif a free, set a=val[b]; else set a free input in b's cone.
        done=False
        for (src,dst) in [(a,b),(b,a)]:
            if dst in freeinp:
                val[dst]=val[src]; changed=True; done=True; break
        if done: continue
        # neither free: set a free input in dst's cone to shift dst by (val[src]-val[dst]) with unit coeff
        for (src,dst) in [(a,b),(b,a)]:
            for w in anc.get(dst,()):
                old=val[w]; val[w]=old+1; forward(); ns['v']=val; d=val[dst]-(val[dst])  # need recompute
                val[w]=old; forward(); ns['v']=val
                break
    forward(); ns['v']=val
    Fq=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
    print(f"iter {it}: {len(lines)-len(Fq)}/{len(lines)} ({len(Fq)} fail)", flush=True)
    if not changed or not Fq: break
forward(); ns['v']=val
Fq=[i for i in range(len(lines)) if eval(eqcode[i],ns)!=0]
print(f"FINAL: {len(lines)-len(Fq)}/{len(lines)} ({len(Fq)} fail): {Fq[:20]}")
if len(Fq)<20:
    json.dump({f"x_{i}":val[i] for i in range(NVARS)}, open('tie_solved.json','w')); print("saved tie_solved.json")
