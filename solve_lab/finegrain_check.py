#!/usr/bin/env python3
"""Verify: are there product-slack atoms where BOTH factors are non-wire (fine-grained)?
A fine-grained slack (product of two non-wire vars, one currently 0 & free) can absorb ANY sub-p value.
Scan all atoms for such slacks and check if any control x_3558 or the gaps."""
import json,pickle
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
mem=defaultdict(int)
for poly in atoms:
    for x in atom_vars(poly): mem[x]+=1
# Find product monomials (a,b) where neither is wire, at least one currently 0
finegrain=[]
for ai,poly in enumerate(atoms):
    for m in poly:
        if len(m)==2:
            a,b=m
            if a not in wire and b not in wire:
                # fine-grained IF one factor is 0 & the other nonzero (activatable to any value)
                va,vb=vA[a],vA[b]
                if (va==0 and vb!=0) or (vb==0 and va!=0):
                    zero = a if va==0 else b
                    nz = b if va==0 else a
                    finegrain.append((ai,zero,nz,vA[nz]%p, zero in freeinp))
print(f"fine-grained activatable product-slacks (non-wire*non-wire, one factor 0): {len(finegrain)}")
# how many have the ZERO factor free (a real knob)?
freeknob=[f for f in finegrain if f[4]]
print(f"  ...with the zero factor being a FREE input: {len(freeknob)}")
# distinct zero-free knobs
zks=set(f[1] for f in freeknob)
print(f"  distinct free zero-knobs: {len(zks)}")
# Do any of these atoms also involve the core (x_3558, x_29322) or gaps (x_4432,x_7068,x_2099,x_19964)?
core_related=[]
targets={3558,29322,4432,7068,2099,19964,642,28730,11150,25739,37758,35389,6671}
for ai,zero,nz,nzval,isfree in freeknob:
    if atom_vars(atoms[ai])&targets:
        core_related.append((ai,zero,nz))
print(f"  fine-grained free-knob slacks in atoms touching core/gaps: {len(core_related)}")
for ai,zero,nz in core_related[:15]:
    print(f"    atom {ai}: activate x_{zero}(free,0) * x_{nz}(={nzval if (nzval:=vA[nz])<10**11 else 'big'}) ; atom touches {sorted(atom_vars(atoms[ai])&targets)}")
