import heal_harness as H, jac_lib as J, time, random, pickle
from collections import defaultdict
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
vd=J.build_duals()
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
Fset=set(F)
NEQ=len(H.eqcode)
# Build Jacobian rows with residuals
rows=[]; consts=[]
for i in range(NEQ):
    rv,grad=J.eq_jac_row(i,vd)
    rows.append(dict(grad)); consts.append(rv)   # consts[i]=residual e_i(base)
# We want certificate c: sum_i c_i * (Jacobian row i) = 0, and sum_i c_i * resid_i != 0.
# Approach: express each FAIL row as combo of SATISFIED rows via sparse elimination WITH provenance.
# Reduce satisfied rows to echelon, provenance = {orig_eq: coef}. Then reduce each fail row; leftover provenance = certificate.
sat=[i for i in range(NEQ) if i not in Fset]
# sparse gaussian on satisfied rows, track provenance
piv_col_row={}     # pivot col -> (rowdict, provdict)
def reduce_row(rowdict, prov):
    rowdict=dict(rowdict); prov=dict(prov)
    while rowdict:
        c=min(rowdict)   # deterministic
        if c in piv_col_row:
            prow,pprov=piv_col_row[c]
            f=rowdict[c]
            for k,v in prow.items():
                nv=(rowdict.get(k,0)-f*v)%p
                if nv: rowdict[k]=nv
                elif k in rowdict: del rowdict[k]
            for k,v in pprov.items():
                nv=(prov.get(k,0)-f*v)%p
                if nv: prov[k]=nv
                elif k in prov: del prov[k]
        else:
            return rowdict,prov,c   # c is a new pivot col
    return rowdict,prov,None
t0=time.time()
for cnt,i in enumerate(sat):
    if not rows[i]:
        continue
    rd,pv,newc=reduce_row(rows[i],{i:1})
    if newc is not None:
        inv=pow(rd[newc],-1,p)
        rd={k:(v*inv)%p for k,v in rd.items()}
        pv={k:(v*inv)%p for k,v in pv.items()}
        piv_col_row[newc]=(rd,pv)
print(f"sat pivots(rank J_sat)={len(piv_col_row)}  time {time.time()-t0:.1f}s",flush=True)
# now reduce each fail row
certs={}
for i in F:
    rd,pv,newc=reduce_row(rows[i],{i:1})
    # if newc is None, fail row fully in sat span; prov pv gives combo: pv (over eqs) has c^T J=0
    if newc is None:
        certs[i]=pv   # this is c with c^T J=0 (rows[i]-combo=0). c^T resid = ?
    else:
        certs[i]=None
        print(f"  fail {i}: row NOT in sat span (independent pivot col {newc})")
# For each cert, compute c^T resid (should be r_i - 0 = resid_i if valid)
print("\nCertificates (fail row expressed in satisfied rowspace):")
usable=[]
for i in F:
    c=certs[i]
    if c is None: continue
    ctr=sum(cf*consts[k] for k,cf in c.items())%p
    print(f"  fail {i}: cert size {len(c)} eqs, c^T resid = {'nonzero' if ctr else 'ZERO'}")
    if ctr: usable.append(i)
pickle.dump({'certs':certs,'consts':consts},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/certs.pkl','wb'))
# CONSERVATION TEST: for a usable cert, Phi(u)=sum c_i e_i(u). Test over random branch-A points.
if usable:
    i0=usable[0]; c=certs[i0]
    print(f"\nConservation test on cert of fail {i0} ({len(c)} eqs):")
    ceqs=list(c.items())
    def Phi():
        ns={'v':H.val,'__builtins__':{}}
        return sum(cf*eval(H.eqcode[k],ns) for k,cf in ceqs)%p
    selectors={4287,2081,9062,20434}
    perturbable=[v for v in sorted(H.freeinp) if v not in selectors]
    for v in H.freeinp: H.val[v]=base[v]
    H.forward(); vals=[Phi()]
    random.seed(21)
    for _ in range(8):
        for v in H.freeinp: H.val[v]=base[v]+(random.randint(-10**9,10**9) if v not in selectors else 0)
        H.forward(); vals.append(Phi())
    for v in H.freeinp: H.val[v]=base[v]
    H.forward()
    allsame=all(x==vals[0] for x in vals)
    print("  Phi CONSERVED (constant)?" , allsame, " base value nonzero?", vals[0]!=0)
    print("  sample values:", [str(x)[:25] for x in vals[:5]])
    if allsame and vals[0]!=0:
        print("  => OBSTRUCTION IS A CONSERVED INVARIANT: branch A PROVABLY INFEASIBLE (2nd-order cannot help)")
    elif not allsame:
        print("  => Phi VARIES: nonlinear escape POSSIBLE, pursue second-order solve")
