#!/usr/bin/env python3
"""Search the 16-ripple equations for both-zero product slacks; attempt second-order heal."""
import heal_harness as H
import json,pickle
from collections import defaultdict
from propagate import load_atoms, atom_vars
p=H.p
atoms=load_atoms()
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
F16=set(H.fails())
print(f"16-ripple state: {len(F16)} fails: {sorted(F16)}")
# For each failing eq, find both-zero free product slacks in it (via atoms mapping is hard; use eqvars + gate structure)
# Simpler: bump PAIRS of free-zero vars in the eq and see if the product term appears (2nd order effect)
import re
VAR=re.compile(r'x_(\d+)')
lines=[L for L in open('../EQUATIONS.txt').read().split('\n') if L.strip()]
mem=defaultdict(int)
for pp in atoms:
    for x in atom_vars(pp): mem[x]+=1
# free-zero vars in the 16 with low membership
freezero=set()
for i in F16:
    for x in H.eqvars[i]:
        if x in H.freeinp and H.val[x]==0 and x not in wire and mem[x]<=8:
            freezero.add(x)
print(f"free-zero low-mem vars in the 16: {len(freezero)}")
# does any PAIR appear as a product in a failing eq (2nd-order slack)? check via the raw eq text for x_a)*(x_b or x_a*x_b patterns
# Use atoms: find atoms with product (a,b) both in freezero, that are in a failing eq's support
so_slacks=[]
fz=freezero
for ai,poly in enumerate(atoms):
    for m in poly:
        if len(m)==2 and m[0] in fz and m[1] in fz and m[0]!=m[1]:
            so_slacks.append((ai,m[0],m[1]))
print(f"both-free-zero product slacks among the 16's free-zero vars: {len(so_slacks)}")
for ai,a,b in so_slacks[:15]:
    print(f"   atom {ai}: x_{a}*x_{b}")
