import heal_harness as H
from propagate import load_atoms, atom_vars
from collections import defaultdict
p=H.p
A=load_atoms()
def evalp(poly,V):
    s=0
    for m,c in poly.items():
        t=c
        for v in m: t*=V[v]
        s+=t
    return s
v013=H.loadd('best/new_instance_partial_39013.json'); vA=H.loadd('best_agentA_39022.json')
changed_free=[2498,2964,4432,6083,7068,11080,14623,14853,23238,24548,28246,31339,36462]
for v in H.freeinp: H.val[v]=v013.get(v,0)
for v in changed_free:
    if v in (4432,7068): continue
    H.val[v]=vA[v]
H.forward()
V=H.val
F=H.fails()
print(f"16-ripple state: {len(F)} fails: {F}")
# nonzero atoms
nz=[(ai,evalp(A[ai],V)) for ai in range(len(A)) if evalp(A[ai],V)!=0]
print(f"nonzero atoms: {len(nz)}: {[ai for ai,_ in nz]}")
# G1, G2 status
print("G1(a20862)=",evalp(A[20862],V), " G2(a20864)=",evalp(A[20864],V))
# a44342, a45677
print("a44342=",evalp(A[44342],V))
print("a45677=",evalp(A[45677],V))
