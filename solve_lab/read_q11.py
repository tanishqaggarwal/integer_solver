#!/usr/bin/env python3
import json,pickle
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
gate_out=set(); gdef={}; consumers=defaultdict(list)
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gate_out.add(dd['t']); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
        for u in dd['vids']: consumers[u].append(dd['t'])
freeinp=set(range(NVARS))-gate_out
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
mem=defaultdict(int)
for pp in atoms:
    for x in atom_vars(pp): mem[x]+=1
d=json.load(open('fc_partial.json'))
v=[0]*NVARS
for k,val in d.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
for ai in [17897,17901,20866,20867,20868,34232]:
    poly=atoms[ai]
    val=sum((lambda m,c:c*__import__('functools').reduce(lambda a,x:a*v[x],m,1))(m,c) for m,c in poly.items())
    print(f"\n=== atom {ai} (val={val if abs(val)<10**11 else 'BIG(%d)'%len(str(abs(val)))}, %p={val%p}) ===")
    for m,c in sorted(poly.items(),key=lambda kv:(len(kv[0]),kv[0])):
        vs="*".join(f"x_{x}" for x in m) if m else "1"
        meta=[f"x{x}[free={x in freeinp},atoms={mem[x]},wire={x in wire},cons={len(consumers[x])},v={v[x] if abs(v[x])<10**9 else 'BIG'}]" for x in m]
        print(f"   {c:+d} * {vs}  {meta}")
