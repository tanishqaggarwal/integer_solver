import os,sys,json,time,random
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import flint
from collections import defaultdict
p=H.p; ctx=flint.fmpz_mod_ctx(p)
pins=json.load(open('pinrec.json'))
allbits=sorted(set(r[1] for r in pins))
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
K-=set(allbits); K=sorted(K); nK=len(K); kidx={v:j for j,v in enumerate(K)}
b=[eval(H.eqcode[i],ns)%p for i in Fp]
# Jacobian m x nK
coldelta=defaultdict(dict); bc={v:H.val[v] for v in K}
for v in K:
    o=bc[v]; setf(v,o+1)
    for r,i in enumerate(Fp):
        d=(eval(H.eqcode[i],ns)-b[r])%p
        if d: coldelta[r][kidx[v]]=d
    setf(v,o)
# Build [J | b] and also J^T augmented with I to get left-null combos
# left null: find c (m-vector) with c^T J =0. Reduce J^T (nK x m) augmented with I_m
JT=flint.fmpz_mod_mat(nK, m+m, ctx)
for r in range(m):
    for c,val in coldelta[r].items(): JT[c,r]=val
for r in range(m): JT[r if r<nK else 0, m+r]=JT[r if r<nK else 0,m+r]  # placeholder
# proper: augment rows differently. Use approach: rref of [J^T | I_m] on nK rows won't track all.
# Simpler: compute rank J and rank[J|b]; get cokernel via nullspace of J^T
Jt=flint.fmpz_mod_mat(nK,m,ctx)
for r in range(m):
    for c,val in coldelta[r].items(): Jt[c,r]=val
rankJ=Jt.rank()
Jm=flint.fmpz_mod_mat(m,nK,ctx)
for r in range(m):
    for c,val in coldelta[r].items(): Jm[r,c]=val
# nullspace of J^T (m-dim vectors c with J^T c =0 i.e. c in left null of J)
NS=Jt.nullspace()  # returns (matrix, nullity) maybe
if isinstance(NS,tuple): Nmat,nul=NS
else: Nmat=NS; nul=Nmat.ncols()
print(f"x_13195=1: {m} fails, rank(J)={rankJ}, cokernel(conserved functionals) dim={m-rankJ}")
# obstruction: for each cokernel vector c, phi = c.b ; count nonzero
bb=flint.fmpz_mod_mat(m,1,ctx)
for r in range(m): bb[r,0]=b[r]
nz_obstr=0; obstr_vecs=[]
for j in range(Nmat.ncols()):
    c=[int(Nmat[r,j]) for r in range(m)]
    phi=sum(c[r]*b[r] for r in range(m))%p
    if phi!=0:
        nz_obstr+=1; obstr_vecs.append((c,phi))
print(f"conserved functionals with NONZERO obstruction (c.b != 0): {nz_obstr}")
# verify conservation of first obstruction vector under random continuous moves
if obstr_vecs:
    c,phi0=obstr_vecs[0]
    supp=[Fp[r] for r in range(m) if c[r]%p!=0]
    print(f"first obstruction: support {len(supp)} eqs: {supp[:15]}")
    def phi():
        return sum(c[r]*(eval(H.eqcode[Fp[r]],ns)) for r in range(m))%p
    vals=[phi()]; random.seed(5); sv={v:H.val[v] for v in K}
    for _ in range(5):
        for v in K: setf(v,sv[v]+random.randint(-3,3))
        vals.append(phi())
        for v in K: setf(v,sv[v])
    print("phi at 6 continuous points:",[str(x)[:12] for x in vals])
    print("CONSERVED?",all(x==vals[0] for x in vals),"nonzero?",vals[0]!=0)
    json.dump({'Fp':Fp,'obstr':[(cc,int(pp)) for cc,pp in obstr_vecs]}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/cok.json','w'))
