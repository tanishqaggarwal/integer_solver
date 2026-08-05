#!/usr/bin/env python3
"""Expanded tangent: does adding core+load continuous DOF resolve the leaf-close obstruction?
Seed duals on a wide column set (leaf-ripple + core + loads + 1-hop closure), gather every
equation with nonzero gradient, solve J.delta = -r mod p, report consistency (inconsistent rows)."""
import heal_harness as H
from jac_lib import D
import flint, sys
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
H.val[17325]=0; H.val[9413]=0; H.forward()
H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]; H.forward()
F1=set(H.fails())
print('failing after leaf-close:',len(F1))

def eqfree(i):
    s=set()
    for v in H.eqvars[i]: s|=H.anc.get(v,{v})
    return s & H.freeinp

# column set: leaf-ripple support + core/load variable ancestors + one-hop
cols=set()
for i in F1: cols|=eqfree(i)
core_load=[29322,3558,33469,27713,35389,6671,11150,25739,37758,14853,12186,24908,16742,1326]
for v in core_load:
    cols|=(H.anc.get(v,{v})&H.freeinp)
# one hop: eqs touching current cols, add their free support (capped)
from collections import defaultdict
var2eq=defaultdict(set)
for i in range(len(H.eqcode)):
    for j in eqfree(i): var2eq[j].add(i)
aff=set()
for j in list(cols): aff|=var2eq[j]
for i in aff: cols|=eqfree(i)
cols=sorted(cols)
print('working columns:',len(cols))
colidx={j:k for k,j in enumerate(cols)}

# seed duals on cols only
val=H.val
vd=[None]*H.NVARS
for j in H.freeinp:
    vd[j]=D(val[j],{colidx[j]:1}) if j in colidx else D(val[j])
ns={'v':vd,'__builtins__':{}}
for k,t in enumerate(H.order):
    r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)

# gather affected equations (nonzero resid or gradient)
rows=[]; rhs=[]; ids=[]
for i in range(len(H.eqcode)):
    rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
    if isinstance(rr,D): rv=rr.v; g=rr.g
    else: rv=rr%p; g={}
    if rv!=0 or g:
        rows.append(g); rhs.append((-rv)%p); ids.append(i)
nr=len(rows); nc=len(cols)
print(f'affected equations {nr}, columns {nc}; failing among them {sum(1 for i in ids if i in F1)}')

# sparse Gaussian elimination mod p (Markowitz-light): eliminate column by column
# Represent rows as dict col->val, plus rhs. Do fraction-free? Use modular inverse.
Rrows=[dict(g) for g in rows]
Rrhs=list(rhs)
pivcol={}   # col -> pivot row index
used=[False]*nr
import sys
def inv(a): return pow(a%p,p-2,p)
order_cols=cols
# choose pivots greedily
for c in range(nc):
    # find a row with nonzero in col c, not used
    pr=-1
    for r in range(nr):
        if not used[r] and Rrows[r].get(c,0)!=0:
            pr=r; break
    if pr<0: continue
    used[pr]=True; pivcol[c]=pr
    iv=inv(Rrows[pr][c])
    # normalize
    Rrows[pr]={k:(v*iv)%p for k,v in Rrows[pr].items()}; Rrhs[pr]=(Rrhs[pr]*iv)%p
    # eliminate from all other rows
    for r in range(nr):
        if r==pr: continue
        f=Rrows[r].get(c,0)
        if f==0: continue
        pv=Rrows[pr]
        rr=Rrows[r]
        for k,v in pv.items():
            nv=(rr.get(k,0)-f*v)%p
            if nv: rr[k]=nv
            elif k in rr: del rr[k]
        Rrhs[r]=(Rrhs[r]-f*Rrhs[pr])%p
# check inconsistency: used rows already reduced; check UNused rows / all rows with empty lhs but nonzero rhs
incons=0; inc_ids=[]
for r in range(nr):
    if not Rrows[r] and Rrhs[r]!=0:
        incons+=1; inc_ids.append(ids[r])
print('rank',len(pivcol),' INCONSISTENT rows',incons)
if incons==0:
    print('*** EXPANDED TANGENT CONSISTENT -> Newton step exists with core/load DOF ***')
else:
    print('still inconsistent; sample inconsistent eq ids:',inc_ids[:10])
