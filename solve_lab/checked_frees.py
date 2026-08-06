#!/usr/bin/env python3
"""Find all equality-check atoms (deg-1 atoms pinning a single free-input leaf) and the frees they pin.
Then: how many free inputs are 'checked' vs genuinely free?"""
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
mem=defaultdict(list)
for ai,poly in enumerate(atoms):
    for v in atom_vars(poly): mem[v].append(ai)
# equality-check atom: deg-1, contains >=1 free input with coeff +-1
checked=set()
checkatom={}
for ai,poly in enumerate(atoms):
    if any(len(m)>1 for m in poly): continue  # deg-1 only
    terms=[(m[0] if m else None,c) for m,c in poly.items()]
    frees_here=[(v,c) for v,c in terms if v is not None and v in freeinp]
    if len(frees_here)==1 and abs(frees_here[0][1])==1:
        v=frees_here[0][0]
        # it's an equality check pinning v IF the other terms are computed (gate outputs)
        checked.add(v); checkatom[v]=ai
print(f"total free inputs: {len(freeinp)}")
print(f"free inputs pinned by an equality-check atom: {len(checked)}")
print(f"genuinely-free (unchecked) inputs: {len(freeinp)-len(checked)}")
# agent A's compensators
compA=[2498,2964,6083,11080,14623,23238,24548,28246,36462,4432,7068]  # non-core-essential changed
print(f"\nagent A compensators checked? {[(v,v in checked) for v in compA]}")
import pickle; pickle.dump({'checked':checked,'checkatom':checkatom,'freeinp':freeinp}, open('checked.pkl','wb'))
