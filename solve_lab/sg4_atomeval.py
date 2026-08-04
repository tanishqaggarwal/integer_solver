import heal_harness as H
import json
from propagate import load_atoms, atom_vars
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
V=H.val
A=load_atoms()
def evalp(poly):
    s=0
    for m,c in poly.items():
        t=c
        for v in m: t*=V[v]
        s+=t
    return s
nz=[]
for ai,poly in enumerate(A):
    r=evalp(poly)
    if r!=0: nz.append((ai,r))
print(f"NONZERO atoms at agentA baseline: {len(nz)}")
for ai,r in nz:
    av=sorted(atom_vars(A[ai]))
    print(f"  atom {ai}: resid={r}  ({r%p} mod p, quot {r//p})  vars={av}")
