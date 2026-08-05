import heal_harness as H, jac_lib as J, time, random, pickle
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
vd=J.build_duals()
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]; Fset=set(F)
NEQ=len(H.eqcode)
rows=[]; consts=[]
for i in range(NEQ):
    rv,grad=J.eq_jac_row(i,vd); rows.append(dict(grad)); consts.append(rv%p)
# Full-reduction Gaussian with provenance. Process ALL eqs; track which give contradictions.
piv={}       # col -> (rowdict, provdict)
def full_reduce(rd, prov):
    rd=dict(rd); prov=dict(prov)
    while True:
        pcol=None
        for c in rd:
            if rd[c]%p and c in piv: pcol=c; break
        if pcol is None: break
        prow,pprov=piv[pcol]; f=rd[pcol]
        for k,v in prow.items():
            nv=(rd.get(k,0)-f*v)%p
            if nv: rd[k]=nv
            elif k in rd: del rd[k]
        for k,v in pprov.items():
            nv=(prov.get(k,0)-f*v)%p
            if nv: prov[k]=nv
            elif k in prov: del prov[k]
    rd={k:v for k,v in rd.items() if v%p}
    return rd, prov
t0=time.time()
certs=[]   # list of (provdict, cval)  where prov^T J = 0 and cval = prov^T resid
# order: satisfied first, then fails, so fails' contradiction provenance is over sat+fail
order=[i for i in range(NEQ) if i not in Fset]+F
for i in order:
    rd,prov=full_reduce(rows[i],{i:1})
    if not rd:
        # row fully reduced to zero in J-part: prov is a left-null vector of J
        cval=sum(cf*consts[k] for k,cf in prov.items())%p
        if cval: certs.append((prov,cval))
    else:
        c=min(rd); inv=pow(rd[c],-1,p)
        piv[c]=({k:(v*inv)%p for k,v in rd.items()},{k:(v*inv)%p for k,v in prov.items()})
print(f"rank(J_all)={len(piv)}  #certificates(left-null with c^T r!=0)={len(certs)}  time {time.time()-t0:.1f}s")
if not certs:
    print("NO first-order obstruction (system first-order CONSISTENT). Newton step exists.")
else:
    # take the sparsest certificate
    certs.sort(key=lambda x:len(x[0]))
    prov,cval=certs[0]
    print(f"sparsest cert: {len(prov)} eqs, c^T r={str(cval)[:20]}..")
    print(f"  eqs in cert:",sorted(prov.keys()))
    # VERIFY c^T J = 0 exactly by direct check
    chk=defaultdict(int)
    for e,cf in prov.items():
        for col,v in rows[e].items(): chk[col]=(chk[col]+cf*v)%p
    bad=sum(1 for v in chk.values() if v%p)
    print(f"  verify c^T J = 0: {'OK' if bad==0 else f'FAIL ({bad} nonzero)'}")
    # CONSERVATION TEST over exact nonlinear branch-A space
    terms=[(e,cf) for e,cf in prov.items() if cf%p]
    def Phi():
        ns={'v':H.val,'__builtins__':{}}
        return sum(cf*eval(H.eqcode[e],ns) for e,cf in terms)%p
    selectors={4287,2081,9062,20434}
    perturbable=[v for v in sorted(H.freeinp) if v not in selectors]
    for v in H.freeinp: H.val[v]=base[v]
    H.forward(); v0=Phi()
    print(f"  Phi(base)={str(v0)[:25]}.. (== c^T r? {v0==cval})")
    # gradient check: perturb each single var by +-1, Phi should be even (grad 0)
    random.seed(7); vals=[v0]; grad_nonzero=0
    for _ in range(6):
        pert={v:random.randint(-10**12,10**12) for v in random.sample(perturbable,150)}
        for v in H.freeinp: H.val[v]=base[v]+pert.get(v,0)
        H.forward(); vp=Phi()
        for v in H.freeinp: H.val[v]=base[v]-pert.get(v,0)
        H.forward(); vm=Phi()
        if (vp-vm)%p!=0: grad_nonzero+=1   # antisymmetric part = 2*grad
        vals.append(vp); vals.append(vm)
    for v in H.freeinp: H.val[v]=base[v]
    H.forward()
    print(f"  gradient nonzero in {grad_nonzero}/6 probes (0 => cert valid, Phi grad=0)")
    print(f"  Phi CONSERVED?", len(set(vals))==1, " distinct values:",len(set(vals)))
    if len(set(vals))==1:
        print("  *** CONSERVED INVARIANT => branch A PROVABLY INFEASIBLE ***")
    else:
        print("  *** Phi VARIES (2nd order) => nonlinear escape possible ***")
    pickle.dump({'terms':terms,'cval':cval},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/phi2.pkl','wb'))
