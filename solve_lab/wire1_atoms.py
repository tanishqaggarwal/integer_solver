#!/usr/bin/env python3
import json
from collections import defaultdict
from propagate import load_atoms, atom_vars, NVARS
p=2**256-2**32-977
atoms=load_atoms()
# signed identity union-find over 2-term |c1|==|c2| atoms (deg-1, 2 vars)
parent,sgn={}, {}
def find(x):
    parent.setdefault(x,x); sgn.setdefault(x,1)
    if parent[x]==x: return x,sgn[x]
    r,s=find(parent[x]); parent[x]=r; sgn[x]=sgn[x]*s; return r,sgn[x]
def union(a,b,rel):
    ra,sa=find(a); rb,sb=find(b)
    if ra!=rb: parent[ra]=rb; sgn[ra]=rel*sb*sa
for poly in atoms:
    if poly.get((),0)!=0: continue
    terms=[(m,c) for m,c in poly.items()]
    if len(terms)==2 and all(len(m)==1 for m,c in terms):
        (m1,c1),(m2,c2)=terms
        if abs(c1)==abs(c2):
            union(m1[0],m2[0], -1 if (c1>0)==(c2>0) else 1)
# class of x_26064 (forced to p)
r26064=find(26064)[0]
wire={y:find(y)[1] for y in list(parent) if find(y)[0]==r26064}
print(f"wire class of x_26064: {len(wire)} members")
for h in [28599,17499,28730,642,26064]:
    print(f"  x_{h} in wire class: {h in wire}  sign={wire.get(h)}")
# load agentA
def loadv(path):
    d=json.load(open(path)); v=[0]*NVARS
    for k,val in d.items():
        idx=int(k[2:]) if k.startswith('x_') else int(k); v[idx]=int(val)
    return v
vA=loadv('best_agentA_39022.json')
# check current wire values (should be p*sign)
wv=set(vA[y]%p for y in wire)
print(f"wire values mod p (distinct): {len(wv)} -> {list(wv)[:3]}")
# SET wire -> sign*1, evaluate atoms
v=vA[:]
for y,s in wire.items(): v[y]=s*1
def ev(poly):
    s=0
    for m,c in poly.items():
        t=c
        for x in m: t*=v[x]
        s+=t
    return s
nz=[ai for ai,poly in enumerate(atoms) if ev(atoms[ai])!=0]
print(f"\nwire->1 (else agentA): NONZERO atoms = {len(nz)}")
# classify: how many involve only wire+partners
pickle_data={'wire':wire,'nz_wire1':nz}
import pickle; pickle.dump(pickle_data, open('wire_data.pkl','wb'))
print("first 40 nonzero atom ids:", nz[:40])
