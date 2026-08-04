#!/usr/bin/env python3
import heal_harness as H
import json,pickle
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=H.p
atoms=load_atoms()
gate_out=set(); gdef={}
with open('atoms/gates.jsonl') as f:
    for line in f:
        dd=json.loads(line); gate_out.add(dd['t']); gdef[dd['t']]=(dd['rhs'],tuple(dd['vids']))
freeinp=set(range(NVARS))-gate_out
mem=defaultdict(int)
for pp in atoms:
    for x in atom_vars(pp): mem[x]+=1
d=json.load(open('sy_regime11_39018.json'))
for v in H.freeinp: H.val[v]=0
for k,vv in d.items(): H.val[int(k[2:])]=int(vv)
H.forward()  # reconstruct
F=set(H.fails())
print(f"sy_39018 reconstructed: {len(F)} fails")
v=H.val
def ev(poly):
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=v[x]
        s+=t
    return s
nz=[ai for ai in range(len(atoms)) if ev(atoms[ai])!=0]
print(f"nonzero atoms: {len(nz)}")
for ai in nz:
    poly=atoms[ai]; deg=max(len(m) for m in poly); val=ev(poly)
    # div-wire: single big coeff * gate_out - product, i.e. coefficient like 13523997 on a var
    bigcoef=[c for m,c in poly.items() if len(m)==1 and abs(c)>10**5]
    freevars=[x for x in atom_vars(poly) if x in freeinp]
    print(f"  atom {ai}: deg{deg} {len(poly)}t val%p={val%p==0 and 'ZEROmodp' or 'nz'} bigcoef={bigcoef[:2]} free-in-atom={freevars[:5]}")
