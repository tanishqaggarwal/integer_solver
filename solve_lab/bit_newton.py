#!/usr/bin/env python3
"""Bit-space move test: is there a direction in the 256 selector bits that reduces the 11
failing residuals while preserving the 39022 satisfied equations (to first order, mod p)?
If consistent -> a lead for discrete search. If inconsistent -> wall confirmed in bit-space too."""
import heal_harness as H
from jac_lib import D
import json, flint
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
F=set(H.fails())
print('failing:',len(F))
selectors=sorted(set(pr[1] for pr in json.load(open('pinrec.json'))) & H.freeinp)
print('selector bits (free):',len(selectors))
bidx={b:k for k,b in enumerate(selectors)}
# seed duals on selector bits only
vd=[None]*H.NVARS
for j in H.freeinp:
    vd[j]=D(H.val[j],{bidx[j]:1}) if j in bidx else D(H.val[j])
ns={'v':vd,'__builtins__':{}}
for k,t in enumerate(H.order):
    r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
# gather rows: all eqs with nonzero grad or residual
rows=[]; rhs=[]; isfail=[]
for i in range(len(H.eqcode)):
    rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
    if isinstance(rr,D): rv=rr.v; g=rr.g
    else: rv=rr%p; g={}
    if rv!=0 or g:
        rows.append(g); rhs.append((-rv)%p); isfail.append(i in F)
nr=len(rows); nc=len(selectors)
print(f'affected-by-bits system: {nr} eqs x {nc} bits; failing among them {sum(isfail)}')
# sparse gaussian elimination mod p to test consistency of J.delta = -r
Rr=[dict(g) for g in rows]; Rb=list(rhs); used=[False]*nr
def inv(a): return pow(a%p,p-2,p)
piv=0
for c in range(nc):
    pr_=-1
    for r in range(nr):
        if not used[r] and Rr[r].get(c,0)!=0: pr_=r;break
    if pr_<0: continue
    used[pr_]=True; piv+=1
    iv=inv(Rr[pr_][c]); Rr[pr_]={k:(v*iv)%p for k,v in Rr[pr_].items()}; Rb[pr_]=(Rb[pr_]*iv)%p
    pv=Rr[pr_]
    for r in range(nr):
        if r==pr_: continue
        f=Rr[r].get(c,0)
        if f==0: continue
        rr_=Rr[r]
        for k,v in pv.items():
            nv=(rr_.get(k,0)-f*v)%p
            if nv: rr_[k]=nv
            elif k in rr_: del rr_[k]
        Rb[r]=(Rb[r]-f*Rb[pr_])%p
incons=[];
for r in range(nr):
    if not Rr[r] and Rb[r]!=0: incons.append(r)
print('rank',piv,' inconsistent rows',len(incons))
if not incons:
    print('*** CONSISTENT: a first-order bit-space direction reduces the 11 fails while preserving the rest ***')
else:
    print('INCONSISTENT in bit-space too (wall confirmed).')
