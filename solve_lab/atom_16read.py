#!/usr/bin/env python3
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
import heal_harness as H
p=2**256-2**32-977
atoms=load_atoms()
mem=defaultdict(int)
for poly in atoms:
    for v in atom_vars(poly): mem[v]+=1
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
D=__import__('pickle').load(open('wire_data.pkl','rb')); wire=D['wire']
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward(); v=H.val
for ai in [7450,7452]:
    print(f"\n=== atom {ai} (GAP) ===")
    for m,c in sorted(atoms[ai].items(),key=lambda kv:(len(kv[0]),kv[0])):
        vs="*".join(f"x_{x}" for x in m) if m else "1"
        meta=[f"x{x}[free={x in freeinp},atoms={mem[x]},wire={x in wire},val={v[x] if abs(v[x])<10**11 else str(v[x])[:6]+'..'}]" for x in m]
        print(f"   {c:+d} * {vs}  {meta}")
# for the deg2 verifiers, show the product terms and the slack partners
for ai in [44342,45677]:
    print(f"\n=== atom {ai} (VERIFIER) product terms ===")
    for m,c in sorted(atoms[ai].items(),key=lambda kv:(len(kv[0]),kv[0])):
        if len(m)==2:
            meta=[f"x{x}[free={x in freeinp},atoms={mem[x]},wire={x in wire},val={v[x] if abs(v[x])<10**11 else 'BIG'}]" for x in m]
            print(f"   {c:+d} * x_{m[0]}*x_{m[1]}  {meta}")
