#!/usr/bin/env python3
"""Key test: can partial_39021's 12 fails be solved using ONLY low-fanout support inputs
(avoiding the cascade)? If the residual is in the column space of low-fanout knobs, solve
cleanly with minimal ripple, iterate Newton, lift, checker."""
import heal_harness as H
from jac_lib import D
import flint, json
p=H.p
d=H.loadd('best/new_instance_partial_39021.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()

def fanout(sinp, exclude):
    # number of eqs (not in exclude) touched by free input sinp
    c=0
    for i in range(len(H.eqcode)):
        if i in exclude: continue
        fa=set()
        for v in H.eqvars[i]: fa|=H.anc.get(v,{v})
        if sinp in fa: c+=1
    return c

def solve_and_step(thresh):
    F=sorted(H.fails())
    supp=set()
    for i in F:
        for v in H.eqvars[i]: supp|=H.anc.get(v,{v})
    supp=(supp & H.freeinp)
    Fset=set(F)
    low=[s for s in supp if fanout(s,Fset)<=thresh]
    # Jacobian of F over low-fanout cols
    colidx={j:k for k,j in enumerate(sorted(H.freeinp))}
    vd=[None]*H.NVARS
    lowset=set(low)
    for j in H.freeinp:
        vd[j]=D(H.val[j],{colidx[j]:1}) if j in lowset else D(H.val[j])
    ns={'v':vd,'__builtins__':{}}
    for k,t in enumerate(H.order):
        r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
    R=[];G=[]
    for i in F:
        rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
        if isinstance(rr,D): R.append(rr.v%p); G.append(rr.g)
        else: R.append(rr%p); G.append({})
    cols=sorted(set().union(*[set(g) for g in G])) if any(G) else []
    sidx={c:k for k,c in enumerate(cols)}
    ctx=flint.fmpz_mod_ctx(p)
    J=flint.fmpz_mod_mat(len(F),len(cols) if cols else 1,ctx)
    JR=flint.fmpz_mod_mat(len(F),(len(cols)+1) if cols else 2,ctx)
    for r,g in enumerate(G):
        for c,co in g.items(): J[r,sidx[c]]=co%p; JR[r,sidx[c]]=co%p
        JR[r,len(cols)]=R[r]%p
    rJ=J.rank(); rJR=JR.rank()
    print(f'  thresh={thresh}: {len(F)} fails, low-fanout cols={len(low)}, active cols={len(cols)}, rank(J)={rJ}, rank[J|R]={rJR}, consistent={rJR==rJ}')
    return F,R,G,cols,sidx,rJ,rJR,colidx

print('=== can low-fanout knobs alone solve the 12? ===')
for th in (0,5,15,30):
    solve_and_step(th)
