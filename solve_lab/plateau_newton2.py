#!/usr/bin/env python3
"""Simultaneous mod-p Newton restricted to the 511-input plateau support (fast).
Seed duals only on the support; gather ALL equations with nonzero residual or gradient;
solve J*d=-R mod p; apply; iterate. Report the affected-eq count (gradient closure) and
whether residuals -> 0 (mod-p solvable) or a conserved obstruction persists."""
import heal_harness as H
from jac_lib import D
import json
p=H.p
PS=json.load(open('plateau_set.json')); supp=PS['supp']
suppset=set(supp); sidx={j:k for k,j in enumerate(supp)}
def inv(a): return pow(a%p,p-2,p)

vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[17325]=0; H.val[9413]=0; H.forward()
H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]
H.val[2964]=H.val[2099]; H.val[24548]=H.val[19964]
for s in (19569,11052):
    if s in H.freeinp: H.val[s]=0
H.forward()
print('aligned fails:',len(H.fails()))

def step():
    vd=[None]*H.NVARS
    for j in H.freeinp:
        vd[j]=D(H.val[j],{sidx[j]:1}) if j in suppset else D(H.val[j])
    ns={'v':vd,'__builtins__':{}}
    for k,t in enumerate(H.order):
        r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
    rows=[];rhs=[];ids=[]
    for i in range(len(H.eqcode)):
        rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
        if isinstance(rr,D): rv=rr.v; g=rr.g
        else: rv=rr%p; g={}
        if rv!=0 or g: rows.append(g); rhs.append((-rv)%p); ids.append(i)
    cols=sorted(set().union(*[set(g) for g in rows])) if rows else []
    Rr=[dict(g) for g in rows]; Rb=list(rhs); used=[False]*len(rows); pivrow={}
    for c in cols:
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
    delta={c:Rb[prr] for c,prr in pivrow.items()}
    return len(rows),len(cols),len(pivrow),incons,delta

for it in range(15):
    nr,nc,rank,incons,delta=step()
    for c,dv in delta.items():
        stp=dv if dv<=p//2 else dv-p
        H.val[supp[c]]=H.val[supp[c]]+stp
    H.forward()
    nf=len(H.fails())
    print(f'it{it}: affected={nr} cols={nc} rank={rank} INCONS={incons} -> fails={nf}', flush=True)
    if nf==0:
        print('*** ZERO FAILS ***'); break
