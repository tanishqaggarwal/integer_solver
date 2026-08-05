#!/usr/bin/env python3
import heal_harness as H
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=H.p
atoms=load_atoms()
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
import pickle
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
vA=H.loadd('best_agentA_39022.json')
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
R12553_old=42775533402728869434716629464193396056515231264222641773817154079369026410240838606908039
gap12553=CONST1-R12553_old
x22972=vA[22972]
print(f"x_12553 fix: x_13458 = x_22972 * x_5081 = gap12553")
print(f"  x_22972 = {x22972}")
print(f"  gap12553 % x_22972 = {gap12553 % x22972} (divisible: {gap12553%x22972==0})")
print(f"  => x_5081 = gap12553 / x_22972 = {gap12553//x22972 if gap12553%x22972==0 else 'NON-INT'}")
# Search ALL atoms containing x_6418 for a non-wire product slack (product where one factor is free&0, other non-wire)
print(f"\n=== x_6418 slack search ===")
for ai,poly in enumerate(atoms):
    if 6418 not in atom_vars(poly): continue
    # products in this atom with a free-zero factor and non-wire other factor
    for m in poly:
        if len(m)==2:
            a,b=m
            for zf,nf in [(a,b),(b,a)]:
                if zf in freeinp and vA[zf]==0 and nf not in wire and nf!=zf:
                    tag="FINE" if abs(vA[nf])>10**20 or vA[nf]==1 else "?"
                    print(f"  atom {ai}: slack x_{zf}(free,0)*x_{nf}(val={vA[nf] if abs(vA[nf])<10**9 else 'BIG'},wire={nf in wire}) [{tag}]")
