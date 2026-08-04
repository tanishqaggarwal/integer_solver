#!/usr/bin/env python3
import json,pickle
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
D=pickle.load(open('wire_data.pkl','rb')); wire=D['wire']
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
v=vA[:]
for y,s in wire.items(): v[y]=s*1
def ev(poly,vv):
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=vv[x]
        s+=t
    return s
# atom membership
from collections import defaultdict
mem=defaultdict(int)
for poly in atoms:
    for vv in atom_vars(poly): mem[vv]+=1
gates=[]
with open('atoms/gates.jsonl') as f:
    for line in f: gates.append(tuple(json.loads(line)['vids']))
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
for ai in [37110,45828]:
    poly=atoms[ai]
    print(f"\n=== atom {ai}: {len(poly)} terms, deg {max(len(m) for m in poly)}, val(wire1)={ev(poly,vv=v)%p} ===")
    for m,c in sorted(poly.items(),key=lambda kv:(len(kv[0]),kv[0])):
        vs="*".join(f"x_{x}" for x in m) if m else "1"
        info=[]
        for x in m:
            info.append(f"x{x}[free={x in freeinp},atoms={mem[x]},wire={x in wire},v={v[x] if v[x]<10**12 else hex(v[x])[:8]+'..'}]")
        print(f"   {c:+d} * {vs}   {info}")
