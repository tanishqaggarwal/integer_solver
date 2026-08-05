#!/usr/bin/env python3
"""Re-audit fullcore_fix with the NEW conserved-functional tool. Old work abandoned it on
fail-count. Measure: is its residual tangent-CONSISTENT (crossable) or does it carry a
conserved obstruction like agentA?"""
import heal_harness as H
from jac_lib import D
import flint
p=H.p

# reproduce fullcore_fix: degenerate core on the 39013 partial
v0=H.loadd('best/new_instance_partial_39013.json')
for v in H.freeinp: H.val[v]=v0.get(v,0)
H.forward()
# set degenerate core (full value)
H.val[14853]=H.val[12186]
H.val[16742]=H.val[24908]
H.forward()
F=set(H.fails())
print('fullcore_fix reproduced: fails =',len(F), sorted(F)[:25])

# check the 2 linear atom conditions
def atomval(terms):
    s=0
    for vl,co in terms:
        t=co
        for x in vl: t=t*H.val[x]
        s+=t
    return s
a25170=H.val[9254]-6788513*(H.val[16742]-H.val[19083])
a27902=H.val[29967]-12846437*(H.val[14853]-H.val[1308])
print('atom25170 =',a25170,' (==0?',a25170==0,')')
print('atom27902 =',a27902,' (==0?',a27902==0,')')
print('freeness: x_9254',9254 in H.freeinp,' x_29967',29967 in H.freeinp,' x_19083',19083 in H.freeinp,' x_1308',1308 in H.freeinp)

# TANGENT CONSISTENCY over ALL free inputs for the failing set
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
    if rv!=0 or g:
        rows.append(g); rhs.append((-rv)%p); ids.append(i)
supp=sorted(set().union(*[set(g) for g in rows]))
sidx={c:k for k,c in enumerate(supp)}
print(f'affected eqs {len(rows)}, support cols {len(supp)}')
# sparse gaussian elimination mod p, consistency
Rr=[dict(g) for g in rows]; Rb=list(rhs); used=[False]*len(rows)
def inv(a): return pow(a%p,p-2,p)
piv=0
for c in supp:
    ci=sidx[c]
    prr=-1
    for r in range(len(rows)):
        if not used[r] and Rr[r].get(c,0)!=0: prr=r;break
    if prr<0: continue
    used[prr]=True; piv+=1
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
incons=[ids[r] for r in range(len(rows)) if not Rr[r] and Rb[r]!=0]
print('rank',piv,' INCONSISTENT rows',len(incons))
if not incons:
    print('*** fullcore_fix residual is TANGENT-CONSISTENT -> a continuous mod-p step closes it! ***')
else:
    print('conserved obstruction remains; sample inconsistent ids:',incons[:8])
