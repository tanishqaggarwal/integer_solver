import heal_harness as H, random
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
def G12():
    return (7376877*H.val[642]+H.val[2099]-H.val[7068])%p, (H.val[4432]-H.val[19964]-H.val[28730])%p
# Test: is G1,G2 exactly linear (affine) in free inputs within branch A?
# affine <=> second finite difference zero: f(x0+a+b)-f(x0+a)-f(x0+b)+f(x0)=0 for all directions.
freelist=sorted(H.freeinp)
random.seed(7)
# keep selectors boolean at their branch-A values by never touching them:
selectors={4287,2081,9062,20434}  # don't perturb
perturbable=[v for v in freelist if v not in selectors]
def setall(assign):
    for v in H.freeinp: H.val[v]=assign[v]
    H.forward()
nonlin_hits=0
for trial in range(20):
    da={v:random.randint(-5,5) for v in random.sample(perturbable,50)}
    db={v:random.randint(-5,5) for v in random.sample(perturbable,50)}
    A=dict(base); 
    for v,x in da.items(): A[v]=base[v]+x
    setall(A); ga=G12()
    B=dict(base)
    for v,x in db.items(): B[v]=base[v]+x
    setall(B); gb=G12()
    AB=dict(base)
    for v,x in da.items(): AB[v]=base[v]+x
    for v,x in db.items(): AB[v]=AB[v]+x
    setall(AB); gab=G12()
    setall(base); g0=G12()
    sd1=(gab[0]-ga[0]-gb[0]+g0[0])%p
    sd2=(gab[1]-ga[1]-gb[1]+g0[1])%p
    if sd1 or sd2: nonlin_hits+=1
print(f"branch-A second-difference nonzero in {nonlin_hits}/20 trials (0 => G1,G2 EXACTLY linear in free inputs)")
# Now certificate conservation in BRANCH A.
# Build local Jacobian over key free vars for the 11 fails, find left-null cert, test Phi conservation.
setall(base)
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
Kset=sorted(set(v for i in F for v in H.eqvars[i] if v in H.freeinp))
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setmany(assign):
    for v,x in assign.items(): H.val[v]=x
    ks=set()
    for v in assign: ks.update(desc_of[v])
    for k in sorted(ks): H.val[H.order[k]]=eval(H.gcode[k],ns)
def residp(): return [eval(H.eqcode[i],ns)%p for i in F]
setmany({v:base[v] for v in H.freeinp})
b0=residp(); inv2=pow(2,p-2,p); m=len(F); nK=len(Kset)
Jm=[[0]*nK for _ in F]
for j,v in enumerate(Kset):
    o=base[v]; setmany({v:o+1}); rp=residp(); setmany({v:o-1}); rm=residp(); setmany({v:o})
    for k in range(m): Jm[k][j]=((rp[k]-rm[k])*inv2)%p
Aug=[Jm[k][:]+[(-b0[k])%p]+[1 if t==k else 0 for t in range(m)] for k in range(m)]
r=0
for c in range(nK+1):
    piv=None
    for i in range(r,m):
        if Aug[i][c]%p: piv=i;break
    if piv is None: continue
    Aug[r],Aug[piv]=Aug[piv],Aug[r]; iv=pow(Aug[r][c]%p,p-2,p)
    Aug[r]=[(x*iv)%p for x in Aug[r]]
    for i in range(m):
        if i!=r and Aug[i][c]%p:
            f=Aug[i][c]; Aug[i]=[(Aug[i][t]-f*Aug[r][t])%p for t in range(len(Aug[i]))]
    r+=1
cert=None
for i in range(m):
    if all(Aug[i][c]%p==0 for c in range(nK)) and Aug[i][nK]%p!=0:
        cert=[Aug[i][nK+1+t]%p for t in range(m)]; break
print("cert found (branch A):",cert is not None)
if cert:
    print("  cert eqs:",[F[k] for k in range(m) if cert[k]%p])
    def Phi():
        rr=[eval(H.eqcode[i],ns) for i in F]
        return sum(cert[k]*rr[k] for k in range(m))%p
    setmany({v:base[v] for v in H.freeinp}); vals=[Phi()]
    random.seed(11)
    for _ in range(8):
        setmany({v:base[v]+(random.randint(-1000,1000) if v not in selectors else 0) for v in H.freeinp})
        vals.append(Phi())
    setmany({v:base[v] for v in H.freeinp})
    print("  Phi over base+8 random branch-A points:", "CONSERVED" if all(x==vals[0] for x in vals) else "VARIES", " nonzero?",vals[0]!=0)
    print("  values:",vals[:4],"...")
