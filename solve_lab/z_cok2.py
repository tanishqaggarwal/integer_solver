import os,sys,json,random
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
from collections import defaultdict
p=H.p
pins=json.load(open('pinrec.json'))
allbits=set(r[1] for r in pins)
vA=H.loadd('best_agentA_39022.json')
base={v:vA.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.val[13195]=1; H.forward()
desc_of=defaultdict(list)
for k,t in enumerate(H.order):
    for w in H.anc[t]: desc_of[w].append(k)
ns={'v':H.val,'__builtins__':{}}
def setf(v,x):
    H.val[v]=x
    for k in desc_of[v]: H.val[H.order[k]]=eval(H.gcode[k],ns)
def mfails(): return [i for i,c in enumerate(H.eqcode) if eval(c,ns)%p!=0]
Fp=mfails(); m=len(Fp)
K=set()
for i in Fp:
    for v in H.eqvars[i]: K|=(H.anc.get(v,{v})&H.freeinp)
K-=allbits; K=sorted(K); nK=len(K)
b=[eval(H.eqcode[i],ns)%p for i in Fp]
bc={v:H.val[v] for v in K}
J=[[0]*nK for _ in range(m)]
for j,v in enumerate(K):
    o=bc[v]; setf(v,o+1)
    for r in range(m): J[r][j]=(eval(H.eqcode[Fp[r]],ns)-b[r])%p
    setf(v,o)
# Augment [ J(nK) | b(1) | I(m) ]; reduce first nK cols; rows with zero J-part & track combo
Aug=[J[r][:]+[b[r]%p]+[1 if t==r else 0 for t in range(m)] for r in range(m)]
def inv(x): return pow(x%p,p-2,p)
piv_r=0
for c in range(nK):
    pr=None
    for r in range(piv_r,m):
        if Aug[r][c]%p: pr=r;break
    if pr is None: continue
    Aug[piv_r],Aug[pr]=Aug[pr],Aug[piv_r]
    iv=inv(Aug[piv_r][c]); Aug[piv_r]=[(x*iv)%p for x in Aug[piv_r]]
    for r in range(m):
        if r!=piv_r and Aug[r][c]%p:
            f=Aug[r][c]; Aug[r]=[(Aug[r][t]-f*Aug[piv_r][t])%p for t in range(len(Aug[r]))]
    piv_r+=1
    if piv_r==m: break
print(f"x_13195=1: {m} fails, rank(J)={piv_r}, cokernel dim={m-piv_r}")
# rows now with zero in first nK cols are cokernel; b-col value = phi (obstruction)
obstr=[]
for r in range(m):
    if all(Aug[r][c]%p==0 for c in range(nK)):
        phi=Aug[r][nK]%p
        c=[Aug[r][nK+1+t]%p for t in range(m)]
        obstr.append((phi,c))
nz=[(phi,c) for phi,c in obstr if phi!=0]
print(f"cokernel functionals: {len(obstr)}, with nonzero obstruction: {len(nz)}")
for phi,c in nz[:5]:
    supp=[Fp[r] for r in range(m) if c[r]%p]
    print(f"  phi={str(phi)[:20]}... support {len(supp)} eqs: {supp}")
# verify conservation of first
if nz:
    phi0,c=nz[0]
    def PHI(): return sum(c[r]*eval(H.eqcode[Fp[r]],ns) for r in range(m))%p
    vals=[PHI()]; random.seed(7); sv={v:H.val[v] for v in K}
    for _ in range(5):
        for v in K: setf(v,sv[v]+random.randint(-4,4))
        vals.append(PHI()); 
        for v in K: setf(v,sv[v])
    print("PHI conserved under continuous moves?",all(x==vals[0] for x in vals),"| nonzero?",vals[0]!=0)
    # obstruction support atoms
    supp=sorted(set(Fp[r] for r in range(m) if c[r]%p))
    print("obstruction eq support:",supp)
