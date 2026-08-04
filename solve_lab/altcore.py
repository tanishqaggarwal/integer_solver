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
v013=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=v013.get(v,0)
H.forward()
r29322 = H.val[29322]%p
r3558 = H.val[3558]%p
print(f"at 39013: x_29322%p={r29322!=0}, x_3558%p={r3558!=0}")
# solve core via x_12186 (for x_29322) and x_16742 (for x_3558)
H.val[12186] = (H.val[12186] + r29322)   # x_29322 = x_14853 - x_12186 -= r29322 => x_12186 += r29322
H.val[16742] = (H.val[16742] - r3558)    # x_3558 = x_24908 - x_16742 += r3558 => x_16742 -= r3558
H.forward()
print(f"after: x_29322%p={H.val[29322]%p!=0}, x_3558%p={H.val[3558]%p!=0}")
F=[i for i in range(len(H.eqcode)) if eval(H.eqcode[i],{'v':H.val,'__builtins__':{}})!=0]
print(f"fails: {len(F)}")
# analyze nonzero atoms
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
pgran=0; finegr=0
for ai in nz:
    poly=atoms[ai]; deg=max(len(m) for m in poly)
    # does it have a fine-grained slack (free-0 * nonwire-nonp)?
    fine=False
    for m in poly:
        if len(m)==2:
            a,b=m
            for zf,nf in [(a,b),(b,a)]:
                if zf in freeinp and v[zf]==0 and nf not in wire and v[nf]!=0 and v[nf]%p!=0:
                    fine=True
    # is it a p-granular gap? deg1 with wire/p slack
    freeres=[x for x in atom_vars(poly) if x in freeinp and abs(v[x])>10**20]
    print(f"  atom {ai}: deg{deg} {len(poly)}t fine-slack={fine} free-res={freeres}")
