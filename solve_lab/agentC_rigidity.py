#!/usr/bin/env python3
"""Decisive test: are x_12186/x_14853 pinned by LINEAR constraints (rigid) or only QUADRATIC product
constraints (2nd-order escapable)? Classify each closure constraint as linear/quadratic in the
handle-moves; compute achievable(da,db,dc) dim using LINEAR-only constraints vs ALL."""
import json, time, ast, re
from agentC_common import (p, gates, order, definer, forward, partial_forward, downstream_ks,
                           val, freeinp, ns, lines, eqcode, eqvars, load_best, CORE, posof, NVARS,
                           rootcode_of, inv)
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
# build closure to ~2500 handles
H=set(CONTROLS); cons=set()
for hop in range(11):
    newcons=set()
    for h in H: newcons|=(affected_eqs(h)-F0)
    cons=newcons
    comp=set()
    for i in cons: comp|=(eqvars[i]&freeinp)
    comp-=deepcone; comp|=set(CONTROLS)
    if comp<=H or len(comp)>2600: break
    H=comp
H=sorted(H); Hidx={h:i for i,h in enumerate(H)}; NH=len(H); cons=sorted(cons)
Hset=set(H)
print(f"handles={NH}, cons={len(cons)}", flush=True)

# classify each constraint eq: quadratic if it has product monomial with BOTH factor-cones
# containing a handle (so the move can create a 2nd-order term); else linear-in-handles.
VAR=re.compile(r'x_(\d+)')
def has_handle_product(expr):
    node=ast.parse(expr,mode='eval').body
    def vids(n): return set(int(m) for m in VAR.findall(ast.unparse(n)))
    found=[False]
    def walk(n):
        if isinstance(n,ast.BinOp):
            if isinstance(n.op,ast.Mult):
                la=vids(n.left); lb=vids(n.right)
                if la and lb:
                    ca=set(); cb=set()
                    for v in la: ca|=freecone(v)
                    for v in lb: cb|=freecone(v)
                    if (ca&Hset) and (cb&Hset):
                        found[0]=True
            walk(n.left); walk(n.right)
        elif isinstance(n,ast.UnaryOp): walk(n.operand)
    walk(node); return found[0]
t0=time.time()
lin_cons=[]; quad_cons=[]
for i in cons:
    if has_handle_product(lines[i].rsplit('=',1)[0]): quad_cons.append(i)
    else: lin_cons.append(i)
print(f"linear-in-handles cons: {len(lin_cons)}; quadratic (handle-product) cons: {len(quad_cons)} ({time.time()-t0:.0f}s)", flush=True)

# Jacobian rows (first-order) for a given constraint subset
def sparse_rank(rows, extra_rows=()):
    pivots={}
    def reduce_row(r):
        r=dict(r)
        while r:
            c=None
            for cc in r:
                if cc in pivots: c=cc; break
            if c is None: return r
            pr=pivots[c]; f=r[c]
            for cc,vv in pr.items():
                r[cc]=(r.get(cc,0)-f*vv)%p
                if r[cc]==0: del r[cc]
        return r
    rank=0
    for r in list(rows)+list(extra_rows):
        rr={c:v for c,v in reduce_row(r).items() if v%p}
        if not rr: continue
        c=next(iter(rr)); ivv=inv(rr[c])
        pivots[c]={cc:(vv*ivv)%p for cc,vv in rr.items()}; rank+=1
    return rank
def build_J(subset):
    base_root={i:eval(rootcode_of(i),ns)%p for i in subset}
    subset_set=set(subset)
    rows_by_eq=defaultdict(dict)
    for h in H:
        dk=downstream_ks(h)
        aff=set(eqbyvar.get(h,()))&subset_set
        for k in dk: aff|=(eqbyvar.get(order[k],())&subset_set)
        o=val[h]; val[h]=o+1; partial_forward(dk)
        for i in aff:
            d=(eval(rootcode_of(i),ns)-base_root[i])%p
            if d: rows_by_eq[i][Hidx[h]]=d
        val[h]=o; partial_forward(dk)
    return [r for r in rows_by_eq.values() if r]
ea={Hidx[12186]:1}; eb={Hidx[14853]:1}; ec={Hidx[16742]:1}
Jlin=build_J(lin_cons)
rL=sparse_rank(Jlin); rLe=sparse_rank(Jlin,[ea,eb,ec])
rLa=sparse_rank(Jlin,[ea]); rLb=sparse_rank(Jlin,[eb]); rLc=sparse_rank(Jlin,[ec])
print(f"\nLINEAR-only constraints: rank={rL}, achievable(da,db,dc) dim={rLe-rL} [a:{rLa-rL} b:{rLb-rL} c:{rLc-rL}]")
Jall=build_J(cons)
rA=sparse_rank(Jall); rAe=sparse_rank(Jall,[ea,eb,ec])
rAa=sparse_rank(Jall,[ea]); rAb=sparse_rank(Jall,[eb]); rAc=sparse_rank(Jall,[ec])
print(f"ALL constraints (1st order): rank={rA}, achievable(da,db,dc) dim={rAe-rA} [a:{rAa-rA} b:{rAb-rA} c:{rAc-rA}]")
print("\nINTERPRETATION:")
print(f"  da movable vs linear-only: {rLa-rL>0}  (True => x_12186 NOT rigidly pinned; only quadratic products pin it => 2nd-order escape exists)")
print(f"  db movable vs linear-only: {rLb-rL>0}")
print(f"  dc movable vs linear-only: {rLc-rL>0}")
