#!/usr/bin/env python3
import heal_harness as H
import json,pickle
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=H.p
atoms=load_atoms()
gate_out=set()
with open('atoms/gates.jsonl') as f:
    for line in f: gate_out.add(json.loads(line)['t'])
freeinp=set(range(NVARS))-gate_out
D=pickle.load(open('wire_data.pkl','rb')); wire=set(D['wire'])
mem=defaultdict(int)
for pp in atoms:
    for x in atom_vars(pp): mem[x]+=1
# load fc_partial
d=json.load(open('fc_partial.json'))
v=[0]*NVARS
for k,val in d.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(val)
def ev(poly):
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=v[x]
        s+=t
    return s
nz=[ai for ai in range(len(atoms)) if ev(atoms[ai])!=0]
print(f"fc_partial (39016, 17 fail): {len(nz)} nonzero atoms")
for ai in nz:
    poly=atoms[ai]; deg=max(len(m) for m in poly)
    # p-granular gap? (deg1 with only wire/p slack) vs fine-grained?
    fine=[]
    pgran=[]
    for m in poly:
        if len(m)==2:
            a,b=m
            for zf,nf in [(a,b),(b,a)]:
                if zf in freeinp and v[zf]==0:
                    if nf in wire or v[nf]%p==0: pgran.append((zf,nf))
                    elif v[nf]!=0: fine.append((zf,nf))
    freeres=[x for x in atom_vars(poly) if x in freeinp and abs(v[x])>10**20]
    print(f"  atom {ai}: deg{deg} {len(poly)}t, FINE-slacks={len(fine)}{fine[:2]}, p-gran={len(pgran)}, free-res={freeres}")
