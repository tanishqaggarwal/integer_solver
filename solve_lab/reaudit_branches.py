#!/usr/bin/env python3
"""Re-audit distinct branches by the NEW metric (conserved-obstruction dimension), not fail-count.
For each saved 39022-class state: close G1/G2 via free leaves, then measure the exact
GF(p) conserved-obstruction dimension on the resulting residual."""
import heal_harness as H
from jac_lib import D
import flint, sys
p=H.p

def measure(path):
    d=H.loadd(path)
    for v in H.freeinp: H.val[v]=d.get(v,0)
    H.forward()
    f0=len(H.fails())
    # close G1/G2 via free leaves (works any branch: sets atoms 20862/20864 to 0)
    H.val[17325]=0; H.val[9413]=0; H.forward()
    H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]; H.forward()
    F=sorted(H.fails())
    if not F:
        print(f'{path}: baseline {f0} -> after close 0 FAILS!!!'); return
    # conserved dimension on F: gradients over ALL free inputs
    freelist=sorted(H.freeinp); colidx={j:k for k,j in enumerate(freelist)}
    vd=[None]*H.NVARS
    for j in H.freeinp: vd[j]=D(H.val[j],{colidx[j]:1})
    ns={'v':vd,'__builtins__':{}}
    for k,t in enumerate(H.order):
        r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
    R=[];G=[]
    for i in F:
        rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
        if isinstance(rr,D): R.append(rr.v%p); G.append(rr.g)
        else: R.append(rr%p); G.append({})
    supp=sorted(set().union(*[set(g) for g in G])) if any(G) else []
    n=len(F); ns_=len(supp); sidx={c:k for k,c in enumerate(supp)}
    ctx=flint.fmpz_mod_ctx(p)
    J=flint.fmpz_mod_mat(n,ns_ if ns_ else 1,ctx)
    JR=flint.fmpz_mod_mat(n,(ns_+1) if ns_ else 2,ctx)
    for r,g in enumerate(G):
        for c,co in g.items(): J[r,sidx[c]]=co%p; JR[r,sidx[c]]=co%p
        JR[r,ns_]=R[r]%p
    rJ=J.rank(); rJR=JR.rank()
    print(f'{path}: baseline {f0} -> after close {n} fails; support {ns_}; rank(J)={rJ} rank[J|R]={rJR} => OBSTRUCTION DIM {rJR-rJ}')

for path in ['best_agentA_39022.json','gadget_handled.json','best/new_instance_partial_39021.json',
             'best/new_instance_partial_39018.json','agentD_route_11.json','agentD_route_01.json',
             'best/new_instance_partial_39007.json']:
    try: measure(path)
    except Exception as e: print(f'{path}: ERROR {type(e).__name__}: {e}')
