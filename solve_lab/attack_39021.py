#!/usr/bin/env python3
"""partial_39021 has obstruction-dim 0 on its residual. Verify full-affected tangent
consistency (ripple included) and attempt an exact mod-p Newton close, then Z-lift + checker."""
import heal_harness as H
from jac_lib import D
import flint
p=H.p
d=H.loadd('best/new_instance_partial_39021.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
H.forward()
F=sorted(H.fails())
print('partial_39021 fails:',len(F),F)

def full_tangent(report=True):
    freelist=sorted(H.freeinp); colidx={j:k for k,j in enumerate(freelist)}
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
    supp=sorted(set().union(*[set(g) for g in rows]))
    if report: print(f'  affected eqs {len(rows)}, support {len(supp)}')
    return rows,rhs,ids,supp,freelist,colidx

rows,rhs,ids,supp,freelist,colidx=full_tangent()
# sparse GE mod p; track solution delta (particular)
def solve(rows,rhs,supp):
    sidx={c:k for k,c in enumerate(supp)}
    Rr=[dict(g) for g in rows]; Rb=list(rhs); used=[False]*len(rows); pivrow={}
    def inv(a): return pow(a%p,p-2,p)
    for c in supp:
        prr=-1
        for r in range(len(rows)):
            if not used[r] and Rr[r].get(c,0)!=0: prr=r;break
        if prr<0: continue
        used[prr]=True; pivrow[c]=prr
        iv=inv(Rr[prr][c]); Rr[prr]={k:(vv*iv)%p for k,vv in Rr[prr].items()}; Rb[prr]=(Rb[prr]*iv)%p
        pv=Rr[prr]
        for r in range(len(rows)):
            if r==prr: continue
            f=Rr[r].get(c,0)
            if f==0: continue
            rr_=Rr[r]
            for k,vv in pv.items():
                nv=(rr_.get(k,0)-f*vv)%p
                if nv: rr_[k]=nv
                elif k in rr_: del rr_[k]
            Rb[r]=(Rb[r]-f*Rb[prr])%p
    incons=[r for r in range(len(rows)) if not Rr[r] and Rb[r]!=0]
    # back out particular solution: free vars=0, pivot var = rb
    delta={}
    for c,prr in pivrow.items():
        delta[c]=Rb[prr]  # since row normalized, pivot var = rhs - sum(nonpiv) ; nonpiv=0
    return incons,delta,pivrow
incons,delta,pivrow=solve(rows,rhs,supp)
print('  INCONSISTENT rows:',len(incons))
if incons:
    print('  not cleanly consistent (ripple obstructs); sample',incons[:5]);
# Newton iterations regardless (nonlinear): apply delta as integer step, forward, recount
import sys
best=len(F)
for it in range(8):
    rows,rhs,ids,supp,freelist,colidx=full_tangent(report=False)
    if not any(r!=0 for r in rhs):
        pass
    incons,delta,pivrow=solve(rows,rhs,supp)
    # apply step: for each pivot column (free-input index), add delta
    for c,dv in delta.items():
        j=freelist[c]
        # center the residue to small signed for stability
        step=dv if dv< p//2 else dv-p
        H.val[j]=(H.val[j]+step)
    H.forward()
    nf=len(H.fails())
    print(f'  Newton it{it}: incons={len(incons)} -> fails={nf}')
    if nf<best: best=nf
    if nf==0:
        print('  *** ZERO FAILS (mod-solve) - checking exact ***'); break
print('best fails reached:',best)
import json
json.dump({str(v):H.val[v] for v in range(H.NVARS)}, open('attack39021_out.json','w'))
print('saved attack39021_out.json')
