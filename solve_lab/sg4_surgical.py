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
def lincoef(ai,v,V):
    s=0
    for m,c in A[ai].items():
        if v in m:
            t=c
            cnt=0
            for w in m:
                if w==v and cnt==0: cnt=1; continue
                t*=V[w]
            s+=t
    return s
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[4287]=1; H.val[9413]=0; H.val[17325]=0
H.forward()
H.val[8731]=H.val[4432]-H.val[20492]
H.val[9118]=H.val[7068]-H.val[37158]
H.forward()
V=H.val
baseF=set(H.fails())
print(f"x_4287=1 + G1/G2 absorbers: {len(baseF)} fails")
# For each currently-nonzero atom, find a fine-grained free absorber (coeff not mult of p, free, currently 0, low fan-out)
def find_absorber(ai):
    r=ev(ai,V)
    if r==0: return None
    best=None
    for v in sorted(atom_vars(A[ai])):
        if v not in H.freeinp: continue
        c=lincoef(ai,v,V)
        if c==0 or c%p==0: continue  # need fine-grained
        if r % c != 0: continue      # need exact integer
        # fan-out: number of OTHER currently-satisfied atoms with v (proxy)
        fo=sum(1 for aj in range(len(A)) if v in atom_vars(A[aj]) and aj!=ai and ev(aj,V)==0)
        if best is None or fo<best[2]:
            best=(v,-r//c, fo, c)  # set val[v] += (-r//c) to zero atom
    return best
import time
t0=time.time()
for rnd in range(6):
    H.forward()
    nz=[ai for ai in range(len(A)) if ev(ai,V)!=0]
    F=set(H.fails())
    print(f"round {rnd}: {len(nz)} nonzero atoms, {len(F)} fails, t={time.time()-t0:.0f}s  atoms={nz[:15]}")
    if not nz: print("ALL ATOMS ZERO"); break
    changed=False
    for ai in nz:
        ab=find_absorber(ai)
        if ab:
            v,delta,fo,c=ab
            V[v]+=delta
            changed=True
    if not changed:
        print("  no fine-grained absorber for remaining atoms:")
        for ai in nz:
            print(f"    a{ai}: r={ev(ai,V)} free-vars-in={[v for v in atom_vars(A[ai]) if v in H.freeinp]}")
        break
H.forward()
F=H.fails()
print(f"FINAL surgical: {len(H.eqcode)-len(F)}/{len(H.eqcode)} ({len(F)} fail)")
