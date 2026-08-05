#!/usr/bin/env python3
"""Solve mod-p with certificate tracking to expose the conservation law pinning atom 41390."""
import pickle, time
p = 2**256-2**32-977
J = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jac2.pkl','rb'))
J_rows=J['J_rows']; constraints=J['constraints']; res=J['res']; ci=J['ci']
def inv(a): return pow(a%p,p-2,p)
# rows with certificate: combo = dict{orig_atom: coef} s.t. row == sum combo*orig_row
rows=[]
for ai in constraints:
    d={c:v%p for c,v in J_rows[ci[ai]].items() if v%p}
    rows.append((d,{ai:1},ai))
pivots={}  # col -> (rowdict without col, combo)
def reduce_row(d,combo):
    d=dict(d); combo=dict(combo)
    stack=[c for c in d if c in pivots]
    while stack:
        col=stack.pop()
        if col not in d: continue
        f=d.pop(col); pd,pcombo=pivots[col]
        for c2,v2 in pd.items():
            nv=(d.get(c2,0)-f*v2)%p
            if nv:
                if c2 not in d and c2 in pivots: stack.append(c2)
                d[c2]=nv
            else: d.pop(c2,None)
        for a2,v2 in pcombo.items():
            nv=(combo.get(a2,0)-f*v2)%p
            if nv: combo[a2]=nv
            else: combo.pop(a2,None)
    return d,combo
incons=[]
for d,combo,ai in rows:
    d,combo=reduce_row(d,combo)
    if not d:
        # check RHS consistency: sum combo*(-res) should be 0
        rhs=sum(cf*(-res[a]) for a,cf in combo.items())%p
        if rhs%p!=0:
            incons.append((ai,combo,rhs))
        continue
    pcol=min(d); ip=inv(d[pcol])
    nd={c:(v*ip)%p for c,v in d.items() if c!=pcol}
    ncombo={a:(v*ip)%p for a,v in combo.items()}
    for c0 in list(pivots):
        pd,pcombo=pivots[c0]
        if pcol in pd:
            f=pd.pop(pcol)
            for c2,v2 in nd.items():
                nv=(pd.get(c2,0)-f*v2)%p
                if nv: pd[c2]=nv
                else: pd.pop(c2,None)
            for a2,v2 in ncombo.items():
                nv=(pcombo.get(a2,0)-f*v2)%p
                if nv: pcombo[a2]=nv
                else: pcombo.pop(a2,None)
            pivots[c0]=(pd,pcombo)
    pivots[pcol]=(nd,ncombo)
print(f"inconsistent: {len(incons)}")
for ai,combo,rhs in incons:
    # combo maps original atoms -> coef; the law is: sum combo*atom_deriv = 0 but sum combo*res != 0
    nzcombo={a:cf for a,cf in combo.items() if a in res and res[a]%p}
    print(f"atom {ai}: forced residual = {rhs}")
    print(f"  conservation law involves gaps (atoms with nonzero res): {nzcombo}")
    print(f"  total constraints in combo: {len(combo)}")
    # verify: rhs == sum combo*res over ALL
    chk=sum(cf*res[a] for a,cf in combo.items())%p
    print(f"  check sum(combo*res)={chk}  (== -rhs mod p? {(chk+rhs)%p==0})")
