#!/usr/bin/env python3
"""At the G1/G2-fixed (x_4432=x_19964, x_7068=x_2099) core-solved state, find nonzero atoms
and analyze their product-slack structure (the wire escape one layer deeper)."""
import json,pickle
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
import heal_harness as H
p=2**256-2**32-977
atoms=load_atoms()
mem=defaultdict(int)
for poly in atoms:
    for v in atom_vars(poly): mem[v]+=1
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
F=H.fails()
print(f"state: {len(F)} eq fails: {F}")
v=H.val
def ev(poly):
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=v[x]
        s+=t
    return s
nz=[(ai,ev(atoms[ai])) for ai in range(len(atoms)) if ev(atoms[ai])!=0]
print(f"nonzero atoms: {len(nz)}")
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
for ai,val in nz:
    poly=atoms[ai]
    vs=sorted(atom_vars(poly))
    deg=max(len(m) for m in poly)
    # product slack terms (deg-2 monomials) with a free/rare factor
    prods=[m for m in poly if len(m)==2]
    rare=[x for x in vs if x in freeinp and mem[x]<=3]
    print(f"atom {ai}: deg{deg}, {len(poly)}t, val%p={val%p}, #prod={len(prods)}, rare-free-partners(<=3 atoms)={rare}")
