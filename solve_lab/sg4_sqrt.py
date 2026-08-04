import heal_harness as H
from propagate import load_atoms, atom_vars
from check_square import try_sqrt
from collections import defaultdict
p=H.p
A=load_atoms()
mem=defaultdict(int)
for poly in A:
    for v in atom_vars(poly): mem[v]+=1
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
for ai in [44342,45677,7450,7452]:
    poly=A[ai]
    deg=max(len(m) for m in poly)
    Q=try_sqrt(poly)
    print(f"=== atom {ai}: deg{deg}, {len(poly)} terms, perfect_square={Q is not None} ===")
    if Q is not None:
        qval=evalp(Q,V)
        print(f"   Q has {len(Q)} terms, Q value = {qval}  (Q^2={qval*qval}, atom={evalp(poly,V)})")
        # linear terms in Q with free & currently-0 vars (fine-grained knobs)
        print("   Q linear terms (deg-1), knobs:")
        for m,c in sorted(Q.items(),key=lambda x:len(x[0])):
            if len(m)==1:
                v=m[0]
                print(f"      {c:+d} * x_{v} [free={v in H.freeinp}, val={V[v] if abs(V[v])<10**11 else 'BIG'}, atoms={mem[v]}]")
