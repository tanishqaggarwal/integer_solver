#!/usr/bin/env python3
"""Track achievable (da,db,dc) dimension as the constraint closure grows to convergence, via
achievable_dim = rank([J; e_a; e_b; e_c]) - rank(J). Decides whether x_12186/x_14853 (da/db) can
move in the wiring null space at all (first order)."""
import json, time, sys
from agentC_common import (p, gates, order, definer, gcode, forward, partial_forward, downstream_ks,
                           val, freeinp, ns, lines, eqcode, eqvars, load_best, CORE, posof, NVARS,
                           pinned, rootcode_of, inv)
from collections import defaultdict

best=load_best(); forward()
gate_defs={t:(rhs,vids) for t,rhs,vids in gates}
_fc={}
def freecone(root):
    if root in _fc: return _fc[root]
    seen=set(); leaves=set(); st=[root]
    while st:
        x=st.pop()
        if x in seen: continue
        seen.add(x)
        if x in gate_defs:
            for u in gate_defs[x][1]: st.append(u)
        elif x in freeinp: leaves.add(x)
    _fc[root]=leaves; return leaves
DEEP=[3558,29322,33469,27713,1326]
deepcone=set()
for d in DEEP+[35389,6671]: deepcone|=freecone(d)
CONTROLS=[12186,14853,16742]
eqbyvar=defaultdict(set)
for i in range(len(lines)):
    for v in eqvars[i]: eqbyvar[v].add(i)
F0=set(i for i in range(len(lines)) if eval(eqcode[i],ns)!=0)
def affected_eqs(h):
    aff=set(eqbyvar.get(h,()))
    for k in downstream_ks(h): aff|=eqbyvar.get(order[k],set())
    return aff

def sparse_rank(rows, extra_rows=()):
    """rows: list of dict{col:val mod p}. Return rank via sparse gaussian elimination mod p."""
    pivots={}   # col -> normalized row dict
    def reduce_row(r):
        r=dict(r)
        while r:
            # leading col = min? any col with pivot -> eliminate
            c=None
            for cc in list(r.keys()):
                if cc in pivots:
                    c=cc; break
            if c is None:
                return r
            pr=pivots[c]; f=r[c]
            for cc,vv in pr.items():
                r[cc]=(r.get(cc,0)-f*vv)%p
                if r[cc]==0: del r[cc]
        return r
    rank=0
    allrows=list(rows)+list(extra_rows)
    for r in allrows:
        rr=reduce_row(r)
        rr={c:v for c,v in rr.items() if v%p}
        if not rr: continue
        c=next(iter(rr)); ivv=inv(rr[c])
        pr={cc:(vv*ivv)%p for cc,vv in rr.items()}
        pivots[c]=pr; rank+=1
    return rank

H=set(CONTROLS); cons=set(); prev=-1
t0=time.time()
for hop in range(12):
    newcons=set()
    for h in H: newcons|=(affected_eqs(h)-F0)
    cons=newcons
    comp=set()
    for i in cons: comp|=(eqvars[i]&freeinp)
    comp-=deepcone; comp|=set(CONTROLS)
    converged = comp<=H
    H2=sorted(comp|H)
    Hidx={h:i for i,h in enumerate(H2)}; NH=len(H2)
    # build Jacobian rows (sparse) for these handles over cons
    dks={h:downstream_ks(h) for h in H2}
    consset=set(cons)
    base_root={i:eval(rootcode_of(i),ns)%p for i in cons}
    Jrows_by_eq=defaultdict(dict)
    for h in H2:
        aff=set(eqbyvar.get(h,()))&consset
        for k in dks[h]: aff|=(eqbyvar.get(order[k],())&consset)
        o=val[h]; val[h]=o+1; partial_forward(dks[h])
        for i in aff:
            d=(eval(rootcode_of(i),ns)-base_root[i])%p
            if d: Jrows_by_eq[i][Hidx[h]]=d
        val[h]=o; partial_forward(dks[h])
    Jrows=[r for r in Jrows_by_eq.values() if r]
    rJ=sparse_rank(Jrows)
    ea={Hidx[12186]:1}; eb={Hidx[14853]:1}; ec={Hidx[16742]:1}
    rJe=sparse_rank(Jrows, extra_rows=[ea,eb,ec])
    achievable=rJe-rJ
    extra=""
    if converged or NH>7000:
        rJa=sparse_rank(Jrows,[ea]); rJb=sparse_rank(Jrows,[eb]); rJc=sparse_rank(Jrows,[ec])
        extra=f" [a:{rJa-rJ} b:{rJb-rJ} c:{rJc-rJ}]"
    print(f"hop {hop}: handles={NH}, cons={len(cons)}, rank(J)={rJ}, achievable(da,db,dc) dim={achievable}"
          f"{extra} converged={converged} ({time.time()-t0:.0f}s)", flush=True)
    if converged: break
    H=comp
    if NH>9000:
        print("handle cap reached; stopping"); break
