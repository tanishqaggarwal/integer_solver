#!/usr/bin/env python3
"""Solve mod-p: drive EVERY sensitive atom to 0 -> b = -res[ai] for all. RREF, test consistency."""
import pickle, time
p = 2**256-2**32-977
J = pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jac2.pkl','rb'))
J_rows=J['J_rows']; constraints=J['constraints']; res=J['res']; ci=J['ci']
def inv(a): return pow(a%p,p-2,p)
rows=[]
for ai in constraints:
    d={c:v%p for c,v in J_rows[ci[ai]].items() if v%p}
    b=(-res[ai])%p
    rows.append((d,b,ai))
pivots={}
def reduce_row(d,b):
    d=dict(d); stack=[c for c in d if c in pivots]
    while stack:
        col=stack.pop()
        if col not in d: continue
        f=d.pop(col); pd,pb=pivots[col]
        for c2,v2 in pd.items():
            nv=(d.get(c2,0)-f*v2)%p
            if nv:
                if c2 not in d and c2 in pivots: stack.append(c2)
                d[c2]=nv
            else: d.pop(c2,None)
        b=(b-f*pb)%p
    return d,b
t0=time.time(); inconsistent=[]
for d,b,ai in rows:
    d,b=reduce_row(d,b)
    if not d:
        if b%p!=0: inconsistent.append((ai,b))
        continue
    pcol=min(d); ip=inv(d[pcol])
    nd={c:(v*ip)%p for c,v in d.items() if c!=pcol}; nb=(b*ip)%p
    for c0 in list(pivots):
        pd,pb=pivots[c0]
        if pcol in pd:
            f=pd.pop(pcol)
            for c2,v2 in nd.items():
                nv=(pd.get(c2,0)-f*v2)%p
                if nv: pd[c2]=nv
                else: pd.pop(c2,None)
            pivots[c0]=(pd,(pb-f*nb)%p)
    pivots[pcol]=(nd,nb)
print(f"elim {time.time()-t0:.1f}s rank={len(pivots)} inconsistent={len(inconsistent)}")
for ai,b in inconsistent[:10]: print(f"   atom {ai}: 0={b}")
if not inconsistent:
    print(">>> CONSISTENT mod p (all sensitive constraints). DOF exists!")
    delta={pcol:pb for pcol,(pd,pb) in pivots.items()}
    nz={c:v for c,v in delta.items() if v%p}
    print(f"particular delta nonzero: {len(nz)}")
    pickle.dump({'delta_modp':delta,'pivots':pivots,'support_free':J['support_free']},
                open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/delta3.pkl','wb'))
    print("saved delta3.pkl")
else:
    pickle.dump({'inconsistent':inconsistent}, open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/incons.pkl','wb'))
