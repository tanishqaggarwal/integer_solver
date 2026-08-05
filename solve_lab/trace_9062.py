#!/usr/bin/env python3
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
gate_out=set(); gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gate_out.add(d['t']); gdef.setdefault(d['t'],(d['rhs'],tuple(d['vids'])))
freeinp=set(range(NVARS))-gate_out
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
# How does x_21279 appear in atoms? alone or in products, and partner values?
print("=== x_21279 appearances in atoms ===")
lone=0; prods=[]
for ai,poly in enumerate(atoms):
    if 21279 not in atom_vars(poly): continue
    for m in poly:
        if 21279 in m:
            if len(m)==1: lone+=1
            elif len(m)==2:
                other=m[0] if m[1]==21279 else m[1]
                prods.append((ai,other,vA[other]))
print(f"  lone (deg-1) appearances: {lone}")
print(f"  product appearances: {len(prods)}; partners & their values:")
pz=defaultdict(int)
for ai,o,ov in prods: pz[o]=ov
for o,ov in pz.items():
    print(f"    x_21279 * x_{o} (partner val={ov if abs(ov)<10**8 else 'BIG'})")
# trace x_9062 upstream: gate + free ancestors
print("\n=== x_9062 gate/upstream ===")
def gate_of(t): return gdef.get(t)
print(f"x_9062 gate: {gate_of(9062)}")
# BFS free ancestors of x_9062
seen=set(); stack=[9062]; fa=set()
while stack:
    t=stack.pop()
    if t in seen: continue
    seen.add(t)
    if t in freeinp: fa.add(t); continue
    if t in gdef:
        for u in gdef[t][1]: stack.append(u)
print(f"x_9062 free ancestors: {len(fa)} -> {sorted(fa)[:15]}")
# Is x_9062 itself a product of things that are 0? show its gate rhs
if 9062 in gdef: print(f"x_9062 rhs: {gdef[9062][0]}")
