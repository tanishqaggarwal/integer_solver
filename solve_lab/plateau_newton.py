#!/usr/bin/env python3
"""The cascade plateaus at ~108 equations -> bounded subsystem. Attack with SIMULTANEOUS
mod-p Newton over the whole plateau free support (not greedy one-at-a-time). Solve J*d=-R
each step over ALL affected equations (so nothing outside breaks). If residuals ->0 mod p,
we have a mod-p solution of the plateau; then lift to Z."""
import heal_harness as H
from jac_lib import D
import flint
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[17325]=0; H.val[9413]=0; H.forward()
H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]
H.val[2964]=H.val[2099]; H.val[24548]=H.val[19964]
for s in (19569,11052):
    if s in H.freeinp: H.val[s]=0
H.forward()
print('aligned fails:',len(H.fails()))

freelist=sorted(H.freeinp); colidx={j:k for k,j in enumerate(freelist)}
def inv(a): return pow(a%p,p-2,p)

def newton_step():
    # seed duals on ALL free inputs; gather affected eqs; solve J d = -R mod p; apply
    vd=[None]*H.NVARS
    for j in H.freeinp: vd[j]=D(H.val[j],{colidx[j]:1})
    ns={'v':vd,'__builtins__':{}}
    for k,t in enumerate(H.order):
        r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
    rows=[];rhs=[];ids=[]
    for i in range(len(H.eqcode)):
        rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
        if isinstance(rr,D): rv=rr.v; g=rr.g
        else: rv=rr%p; g={}
        if rv!=0 or g: rows.append(g); rhs.append((-rv)%p); ids.append(i)
    supp=sorted(set().union(*[set(g) for g in rows])) if rows else []
    # sparse GE with pivot tracking for particular solution (min-support)
    Rr=[dict(g) for g in rows]; Rb=list(rhs); used=[False]*len(rows); pivrow={}
    for c in supp:
        prr=-1
        for r in range(len(rows)):
            if not used[r] and Rr[r].get(c,0)!=0: prr=r;break
        if prr<0: continue
        used[prr]=True; pivrow[c]=prr
        iv=inv(Rr[prr][c]); Rr[prr]={k:(v*iv)%p for k,v in Rr[prr].items()}; Rb[prr]=(Rb[prr]*iv)%p
        pv=Rr[prr]
        for r in range(len(rows)):
            if r==prr: continue
            f=Rr[r].get(c,0)
            if f==0: continue
            rr_=Rr[r]
            for k,v in pv.items():
                nv=(rr_.get(k,0)-f*v)%p
                if nv: rr_[k]=nv
                elif k in rr_: del rr_[k]
            Rb[r]=(Rb[r]-f*Rb[prr])%p
    incons=sum(1 for r in range(len(rows)) if not Rr[r] and Rb[r]!=0)
    # particular solution: pivot var = its reduced rhs (free non-pivot vars stay 0)
    delta={}
    for c,prr in pivrow.items(): delta[c]=Rb[prr]
    return len(rows),len(supp),len(pivrow),incons,delta

for it in range(12):
    nr,nc,rank,incons,delta=newton_step()
    # apply delta (centered) as integer step
    for c,dv in delta.items():
        step=dv if dv<=p//2 else dv-p
        H.val[freelist[c]]=H.val[freelist[c]]+step
    H.forward()
    nf=len(H.fails())
    print(f'it{it}: affected={nr} supp={nc} rank={rank} incons={incons} -> fails={nf}')
    if nf==0:
        print('*** ZERO FAILS mod-solve; verifying exact ***')
        break
