import heal_harness as H, jac_lib as J, time, random
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
vd=J.build_duals()
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
Fset=set(F); NEQ=len(H.eqcode)
rows=[]; consts=[]
for i in range(NEQ):
    rv,grad=J.eq_jac_row(i,vd); rows.append(dict(grad)); consts.append(rv%p)
sat=[i for i in range(NEQ) if i not in Fset]
# sat echelon WITH provenance
piv={}  # col -> (rowdict, provdict{eq:coef})
def reduce_full(rowdict,prov):
    rowdict=dict(rowdict); prov=dict(prov)
    while rowdict:
        c=min(rowdict)
        if c in piv:
            prow,pprov=piv[c]; f=rowdict[c]
            for k,v in prow.items():
                nv=(rowdict.get(k,0)-f*v)%p
                if nv: rowdict[k]=nv
                elif k in rowdict: del rowdict[k]
            for k,v in pprov.items():
                nv=(prov.get(k,0)-f*v)%p
                if nv: prov[k]=nv
                elif k in prov: del prov[k]
        else:
            return rowdict,prov,c
    return rowdict,prov,None
t0=time.time()
for i in sat:
    if not rows[i]: continue
    rd,pv,nc=reduce_full(rows[i],{i:1})
    if nc is not None:
        inv=pow(rd[nc],-1,p)
        piv[nc]={k:(v*inv)%p for k,v in rd.items()},{k:(v*inv)%p for k,v in pv.items()}
        piv[nc]=(piv[nc][0],piv[nc][1])
print(f"sat rank={len(piv)} time {time.time()-t0:.1f}s")
# reduce each fail row: get residual R_i (leftover cols) and provenance P_i (sat combo)
R=[]; P=[]
for i in F:
    rd,pv,nc=reduce_full(rows[i],{})   # prov starts empty: track only sat coefs used
    R.append(rd); P.append(pv)
# find lambda (11-vec) with sum lambda_i R_i = 0.  Stack R_i as rows, find left-null.
allcols=sorted(set().union(*[set(r) for r in R]))
print("leftover cols after full sat-reduction:",len(allcols),"->",allcols[:8],"...")
M=[[R[i].get(c,0) for c in allcols] for i in range(len(F))]  # 11 x len(allcols)
# left null space: lambda with lambda^T M = 0. Do it by augmenting with identity and reducing columns.
m=len(F); ncol=len(allcols)
A=[M[i][:]+[1 if t==i else 0 for t in range(m)] for i in range(m)]
r=0; 
for c in range(ncol):
    piv2=None
    for ii in range(r,m):
        if A[ii][c]%p: piv2=ii;break
    if piv2 is None: continue
    A[r],A[piv2]=A[piv2],A[r]; iv=pow(A[r][c]%p,-1,p)
    A[r]=[(x*iv)%p for x in A[r]]
    for ii in range(m):
        if ii!=r and A[ii][c]%p:
            f=A[ii][c]; A[ii]=[(A[ii][t]-f*A[r][t])%p for t in range(len(A[ii]))]
    r+=1
print(f"rank of R (fails mod sat) = {r}, so left-null dim = {m-r}")
lams=[]
for ii in range(r,m):
    if all(A[ii][c]%p==0 for c in range(ncol)):
        lam=[A[ii][ncol+t]%p for t in range(m)]
        lams.append(lam)
print(f"found {len(lams)} certificate lambdas")
# pick lambda with sum lam_i * consts[F[i]] != 0
chosen=None
for lam in lams:
    ctr=sum(lam[i]*consts[F[i]] for i in range(m))%p
    if ctr: chosen=(lam,ctr);break
if chosen is None:
    print("all certificate lambdas give c^T r = 0 -> residual IS in first-order-reachable set?!")
else:
    lam,ctr=chosen
    print(f"chosen lambda, c^T r = {str(ctr)[:20]}.. (nonzero => genuine first-order obstruction)")
    # build full certificate provenance: c = sum lam_i (e_{F[i]} - P_i . sat)
    mu=defaultdict(int)
    for i in range(m):
        for e,cf in P[i].items(): mu[e]=(mu[e]+lam[i]*cf)%p
    # Phi(u) = sum lam_i e_{F[i]}(u) - sum mu_e e_e(u)
    terms=[(F[i],lam[i]) for i in range(m) if lam[i]]+[(e,(-cf)%p) for e,cf in mu.items() if cf]
    print(f"certificate Phi uses {len(terms)} equations")
    def Phi():
        ns={'v':H.val,'__builtins__':{}}
        return sum(cf*eval(H.eqcode[e],ns) for e,cf in terms)%p
    # sanity: gradient zero at base -> Phi second-diff should reveal conservation
    selectors={4287,2081,9062,20434}
    perturbable=[v for v in sorted(H.freeinp) if v not in selectors]
    for v in H.freeinp: H.val[v]=base[v]
    H.forward(); v0=Phi()
    print("Phi(base) =",str(v0)[:30],"(should equal c^T r =",str(ctr)[:15],")")
    random.seed(99); vals=[v0]; firstorder_ok=True
    # test first-order: tiny perturbation should NOT change Phi (grad 0)
    for scale,label in [(1,'tiny'),(1,'tiny'),(10**6,'mid'),(10**18,'big'),(10**40,'huge')]:
        for _ in range(2):
            for v in H.freeinp: H.val[v]=base[v]+(random.randint(-scale,scale) if v not in selectors else 0)
            H.forward(); vals.append(Phi())
    for v in H.freeinp: H.val[v]=base[v]
    H.forward()
    allsame=all(x==v0 for x in vals)
    print("Phi CONSERVED over all probes?",allsame)
    print("distinct Phi values:",len(set(vals)))
    for x in vals[:8]: print("   ",str(x)[:40])
    if allsame:
        print("\n*** OBSTRUCTION IS A CONSERVED INVARIANT (Phi=const!=0) => branch A PROVABLY INFEASIBLE; no 2nd/higher-order move helps ***")
    else:
        print("\n*** Phi VARIES => the invariant is breakable; nonlinear escape POSSIBLE ***")
        import pickle
        pickle.dump({'terms':terms,'ctr':ctr},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/phi.pkl','wb'))
