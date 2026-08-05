import heal_harness as H
from propagate import load_atoms, atom_vars
from collections import defaultdict
p=H.p
A=load_atoms()
mem=defaultdict(int)
for poly in A:
    for v in atom_vars(poly): mem[v]+=1
def ev(ai,V):
    s=0
    for m,c in A[ai].items():
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
F=set(H.fails())
# which failing eq contains which residual atom - map atom->eqs via eqs it appears in
# Actually map: for each failing eq, which of the 4 atoms explain it
resid_atoms=[7450,7452,44342,45677]
print("residual atom values:")
for ai in resid_atoms:
    print(f"  a{ai} = {ev(ai,V)}   (mod p = {ev(ai,V)%p})")
# a7450, a7452 healers: free ancestors of their gate vars
for ai in [7450,7452]:
    print(f"\n=== a{ai} healers ===")
    for v in sorted(atom_vars(A[ai])):
        anc=H.anc.get(v,set())
        fa=sorted(a for a in anc if a in H.freeinp)
        print(f"  x_{v}: free={v in H.freeinp} val={V[v] if abs(V[v])<1e12 else 'BIG'} #freeanc={len(fa)} {'freeanc='+str(fa[:8]) if not v in H.freeinp else ''}")
# a44342, a45677: extract linear coefficients wrt candidate fine-grained knobs
print("\n=== a44342, a45677 linear structure (coeff of each knob) ===")
def lincoef(ai, v):
    # d(atom)/d(v) at current point = sum over monomials containing v of c*prod(other factors)
    s=0
    for m,c in A[ai].items():
        if v in m:
            t=c*m.count(v)  # (v appears once typically)
            for w in m:
                if w==v: continue
                t*=V[w]
            s+=t
    return s
knobs=[9280,27711,30175,16900,34661,36215,21488,1786,6497,18312,6831,11052,19569,29062,19126,22511]
for ai in [44342,45677]:
    print(f"a{ai} (resid={ev(ai,V)%p} mod p):")
    for v in knobs:
        c=lincoef(ai,v)
        if c!=0:
            print(f"   d/dx_{v} = {c if abs(c)<1e12 else str(c%p)+' (modp)'}  [free={v in H.freeinp}, val={V[v]}, atoms={mem[v]}]")
