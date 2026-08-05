#!/usr/bin/env python3
import heal_harness as H
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=H.p
atoms=load_atoms()
gate_out=set(); gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        d=json.loads(line); gate_out.add(d['t']); gdef[d['t']]=(d['rhs'],tuple(d['vids']))
freeinp=set(range(NVARS))-gate_out
mem=defaultdict(int)
for pp in atoms:
    for x in atom_vars(pp): mem[x]+=1
vA=H.loadd('best_agentA_39022.json')
CONST1=97171863764434070215824145711260403004952728652948669662983319257693684265837195009100680
CONST2=126767545623909574255290391153759363968073470399639361054829680359428658595949132261910506
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[12553]=CONST1; H.val[6418]=CONST2
H.forward()
v=H.val
def ev(poly):
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=v[x]
        s+=t
    return s
nz=[ai for ai in range(len(atoms)) if ev(atoms[ai])!=0]
print(f"nonzero atoms at new state: {len(nz)}")
for ai in nz:
    poly=atoms[ai]; vs=sorted(atom_vars(poly)); deg=max(len(m) for m in poly)
    # find free residue in it (like x_12553): a free input whose value is BIG and coeff +-1
    freeres=[x for x in vs if x in freeinp and abs(v[x])>10**20]
    print(f"  atom {ai}: deg{deg} {len(poly)}t, free-BIG-residues={freeres}")
