#!/usr/bin/env python3
import json
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
v=loadv('best_agentA_39022.json')
v013=loadv('best/new_instance_partial_39013.json')
for ai in [20862,20864,42669]:
    poly=atoms[ai]
    print(f"\n=== atom {ai} ===")
    for m,c in sorted(poly.items(), key=lambda kv:(len(kv[0]),kv[0])):
        vs="*".join(f"x_{x}" for x in m) if m else "1"
        print(f"   {c:+d} * {vs}    [vals agentA: {[v[x] for x in m]}]")
    # value
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=v[x]
        s+=t
    print(f"   value at agentA = {s}  (%p={s%p})")
    # value at 39013
    s2=0
    for m,c in poly.items():
        t=c
        for x in m: t*=v013[x]
        s2+=t
    print(f"   value at 39013  = {s2}  (%p={s2%p})")
