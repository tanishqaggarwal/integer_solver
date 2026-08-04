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
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
mem=defaultdict(int)
for pp in atoms:
    for x in atom_vars(pp): mem[x]+=1
d=json.load(open('fc_partial.json'))
for v in H.freeinp: H.val[v]=d.get('x_%d'%v,0) if isinstance(d.get('x_%d'%v,0),int) else int(d.get('x_%d'%v,0))
for k,vv in d.items(): H.val[int(k[2:])]=int(vv)
# apply gadget fix
H.val[14865]=H.val[12553]; H.val[31861]=H.val[6418]; H.val[33168]=0
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
print(f"29-fail state: {len(nz)} nonzero atoms")
for ai in nz:
    poly=atoms[ai]; deg=max(len(m) for m in poly)
    # p-granular? (wire/p factor slack) vs fine?
    pg=any((len(m)==2 and (m[0] in wire or m[1] in wire or v[m[0]]%p==0 or v[m[1]]%p==0)) for m in poly if len(m)==2)
    # private free slack?
    priv=[x for x in atom_vars(poly) if x in freeinp and mem[x]<=4]
    print(f"  atom {ai}: deg{deg} {len(poly)}t, has-wire/p-factor={pg}, private-free={priv[:4]}")
