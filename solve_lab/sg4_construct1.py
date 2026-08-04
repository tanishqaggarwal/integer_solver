import heal_harness as H
from propagate import load_atoms, atom_vars
from collections import defaultdict
p=H.p
A=load_atoms()
mem=defaultdict(int)
for poly in A:
    for v in atom_vars(poly): mem[v]+=1
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
baseF=set(H.fails())
V0=H.val[:]
# STEP: activate x_4287=1, forward to get shifted bases
H.val[4287]=1
H.val[9413]=0; H.val[17325]=0
H.forward()
x20492=H.val[20492]; x37158=H.val[37158]
x4432=H.val[4432]; x7068=H.val[7068]
# gaps post-flip
g2=x4432-x20492
g1=x7068-x37158
print("post-flip x_20492=",str(x20492)[:12],"  x_37158=",str(x37158)[:12])
print("g2(=x_8731)=",g2)
print("g1(=x_9118)=",g1)
# set the fine-grained absorbers
H.val[8731]=g2
H.val[9118]=g1
H.forward()
# verify G1,G2
print("\nG1(a20862)=",sum(c*__import__('functools').reduce(lambda a,b:a*H.val[b],m,1) for m,c in A[20862].items()))
def ev(ai):
    s=0
    for m,c in A[ai].items():
        t=c
        for v in m: t*=H.val[v]
        s+=t
    return s
print("G1=",ev(20862)," G2=",ev(20864)," verifier=",ev(42669))
F=set(H.fails())
print(f"\nfails: {len(F)}  (baseline was {len(baseF)})")
newbroke=sorted(F-baseF); fixed=sorted(baseF-F)
print("newly broke:",len(newbroke), newbroke[:40])
print("fixed:",len(fixed),fixed)
# nonzero atoms now
nz=[ai for ai in range(len(A)) if ev(ai)!=0]
print(f"nonzero atoms: {len(nz)}: {nz[:40]}")
