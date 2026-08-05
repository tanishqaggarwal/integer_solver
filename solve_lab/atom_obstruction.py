#!/usr/bin/env python3
"""METHODOLOGY method: evaluate all atoms at agentA's 39022 solution, find nonzero (obstruction) atoms,
read them raw to find product-slack absorbers."""
import json,time
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
t0=time.time()
atoms=load_atoms()
print(f"loaded {len(atoms)} atoms in {time.time()-t0:.1f}s")
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
v=loadv('best_agentA_39022.json')
def ev(poly):
    s=0
    for m,c in poly.items():
        t=c
        for var in m: t*=v[var]
        s+=t
    return s
nz=[]
for ai,poly in enumerate(atoms):
    val=ev(poly)
    if val!=0: nz.append((ai,val))
print(f"NONZERO atoms at agentA 39022: {len(nz)}")
for ai,val in nz:
    poly=atoms[ai]
    vs=sorted(atom_vars(poly))
    deg=max(len(m) for m in poly)
    print(f"  atom {ai}: deg={deg}, {len(poly)} terms, {len(vs)} vars, val%p={val%p}, val==0mod p:{val%p==0}")
import pickle
pickle.dump({'nz':nz}, open('atom_nz.pkl','wb'))
