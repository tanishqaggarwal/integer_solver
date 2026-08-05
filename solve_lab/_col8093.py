import heal_harness as H, jac_lib as J, pickle
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
print("free var at column 8093 =",J.freelist[8093])
vd=J.build_duals()
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]
# reduce each fail row against sat rowspace, get (c_i coefficient on col8093, residual r_i), check proportionality
D=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/certs.pkl','rb'))
consts=D['consts']
# rebuild sat pivot structure to reduce fail rows and read the col-8093 coefficient
from collections import defaultdict
NEQ=len(H.eqcode); Fset=set(F)
rows=[]
for i in range(NEQ):
    rv,grad=J.eq_jac_row(i,vd); rows.append(dict(grad))
sat=[i for i in range(NEQ) if i not in Fset]
piv={}
def reduce_row(rowdict):
    rowdict=dict(rowdict)
    while rowdict:
        c=min(rowdict)
        if c in piv:
            prow=piv[c]; f=rowdict[c]
            for k,v in prow.items():
                nv=(rowdict.get(k,0)-f*v)%p
                if nv: rowdict[k]=nv
                elif k in rowdict: del rowdict[k]
        else:
            return rowdict,c
    return rowdict,None
for i in sat:
    if not rows[i]: continue
    rd,nc=reduce_row(rows[i])
    if nc is not None:
        inv=pow(rd[nc],-1,p); rd={k:(v*inv)%p for k,v in rd.items()}; piv[nc]=rd
print("\nfail: residual r_i, coeff on col8093 after sat-reduction, ratio -r_i/c_i:")
ratios=[]
for i in F:
    rd,nc=reduce_row(rows[i])
    c8093=rd.get(8093,0)
    ri=consts[i]%p
    if c8093:
        ratio=(-ri)*pow(c8093,-1,p)%p
        ratios.append(ratio)
    else:
        ratio=None
    # also show other leftover cols
    leftover={k:v for k,v in rd.items()}
    print(f"  fail {i}: r={'0' if ri==0 else str(ri)[:12]+'..'}, leftover cols={sorted(leftover)}, ratio={'--' if ratio is None else str(ratio)[:20]+'..'}")
uniq=set(r for r in ratios if r is not None)
print(f"\ndistinct ratios among fails with col8093: {len(uniq)}")
print("=> FIRST-ORDER CONSISTENT along col8093?" , len(uniq)<=1)
