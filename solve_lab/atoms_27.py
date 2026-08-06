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
vA=H.loadd('best_agentA_39022.json'); v013=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=v013.get(v,0)
H.val[14853]=vA[14853]; H.val[31339]=vA[31339]
H.forward(); v=H.val
def ev(poly):
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=v[x]
        s+=t
    return s
nz=[ai for ai in range(len(atoms)) if ev(atoms[ai])!=0]
print(f"core-only state (27 fails): {len(nz)} nonzero atoms")
for ai in nz:
    poly=atoms[ai]; vs=sorted(atom_vars(poly)); deg=max(len(m) for m in poly)
    # classify: gap (deg1) or verifier (deg2+); free-residues; fine-grained slack (both-0 free non-wire product)
    freeres=[x for x in vs if x in freeinp and abs(v[x])>10**20]
    finegrain=[]
    for m in poly:
        if len(m)==2:
            a,b=m
            for zf,nf in [(a,b),(b,a)]:
                if zf in freeinp and v[zf]==0 and nf not in wire and (v[nf]==1 or (v[nf]!=0 and v[nf]%p!=0)):
                    finegrain.append((zf,nf))
    print(f"  atom {ai}: deg{deg} {len(poly)}t, free-residues={freeres}, fine-slacks={finegrain[:3]}")
