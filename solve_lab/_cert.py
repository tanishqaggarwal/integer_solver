import heal_harness as H, re, random
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
base[4287]=1; base[2081]=1; base[9118]=base[7068]; base[8731]=base[4432]
for v in H.freeinp: H.val[v]=base[v]
H.forward()
F=H.fails()
Kset=sorted(set(v for i in F for v in H.eqvars[i] if v in H.freeinp))
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
touch={v:sorted(set(desc_of[v])) for v in Kset}
def setmany(assign):
    for v,x in assign.items(): H.val[v]=x
    ks=set()
    for v in assign: ks.update(touch[v])
    for k in sorted(ks): H.val[H.order[k]]=eval(H.gcode[k],ns)
def residp():
    return [eval(H.eqcode[i],ns)%p for i in F]
b0=residp(); inv2=pow(2,p-2,p); m=len(F); nK=len(Kset)
J=[[0]*nK for _ in F]
for j,v in enumerate(Kset):
    o=base[v]; setmany({v:o+1}); rp=residp(); setmany({v:o-1}); rm=residp(); setmany({v:o})
    for k in range(m): J[k][j]=((rp[k]-rm[k])*inv2)%p
# find left-null certificate: c such that c^T J = 0 but c^T b != 0
# augment [J | -b0] and row-reduce tracking the combination
Aug=[J[k][:]+[(-b0[k])%p]+[1 if t==k else 0 for t in range(m)] for k in range(m)]  # last m cols track c
# reduce first nK columns; if a row has all-zero in first nK+1 but nonzero... find inconsistency
cols=nK+1
r=0; pivrows=[]
for c in range(cols):
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
# now find a row with zeros in first nK cols but nonzero in col nK (the -b) => inconsistent; its c is cols nK+1..end
cert=None
for i in range(m):
    if all(Aug[i][c]%p==0 for c in range(nK)) and Aug[i][nK]%p!=0:
        cert=[Aug[i][nK+1+t]%p for t in range(m)]; break
print("certificate found:", cert is not None)
if cert:
    nzc=[(F[k],cert[k]) for k in range(m) if cert[k]%p]
    print(f"certificate uses {len(nzc)} equations: {[F[k] for k in range(m) if cert[k]%p]}")
    # Phi = sum c_k * resid_k  (mod p). Test conservation at random points.
    def Phi():
        rr=[eval(H.eqcode[i],ns) for i in F]
        return sum(cert[k]*rr[k] for k in range(m))%p
    setmany({v:base[v] for v in Kset})
    vals=[Phi()]
    random.seed(3)
    for _ in range(6):
        setmany({v:base[v]+random.randint(-1000,1000) for v in Kset})
        vals.append(Phi())
    setmany({v:base[v] for v in Kset})
    print("Phi at base + 6 random points (mod p):")
    for x in vals: print("  ",x)
    print("CONSERVED (all equal)?", all(x==vals[0] for x in vals), " nonzero?", vals[0]!=0)
