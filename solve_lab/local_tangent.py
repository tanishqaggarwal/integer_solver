#!/usr/bin/env python3
"""Local tangent solve: after setting leaves to CONSTs (closing G1/G2), find a mod-p
step in the leaf-ripple free inputs that zeros all affected equations. Restrict duals
to the leaf-ripple columns so gradients are cheap."""
import heal_harness as H
import jac_lib
from jac_lib import D
import flint
p=H.p

# 1. baseline agentA
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
F0=set(H.fails())
print('baseline fails',len(F0))

# 2. close G1/G2 by setting leaves = computed x_2099/x_19964, slacks 0
H.val[17325]=0; H.val[9413]=0
H.forward()
C2099=H.val[2099]; C19964=H.val[19964]
H.val[7068]=C2099; H.val[4432]=C19964
H.forward()
F1=set(H.fails())
print('after leaf-close fails',len(F1),'  fixed',len(F0-F1),' broke',sorted(F1-F0))

# 3. working columns = free inputs in leaf-ripple support (+ leaves/slacks)
leaf=sorted(F1)  # the current failing equations
supp=set()
for i in leaf:
    for v in H.eqvars[i]: supp|=H.anc.get(v,{v})
cols=sorted((supp & H.freeinp))
print('working free columns:',len(cols))
colidx={j:k for k,j in enumerate(cols)}

# 4. build duals seeded ONLY on cols
val=H.val
vd=[None]*H.NVARS
for j in H.freeinp:
    if j in colidx: vd[j]=D(val[j],{colidx[j]:1})
    else: vd[j]=D(val[j])
ns={'v':vd,'__builtins__':{}}
for k,t in enumerate(H.order):
    r=eval(H.gcode[k],ns)
    vd[t]=r if isinstance(r,D) else D(r)

# 5. for ALL equations, get residual + grad over cols; keep rows that are affected
rows=[]; rhs=[]; affected=[]
for i in range(len(H.eqcode)):
    rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
    if isinstance(rr,D):
        rv=rr.v; g=rr.g
    else:
        rv=rr%p; g={}
    if rv!=0 or g:
        rows.append(g); rhs.append((-rv)%p); affected.append(i)
print('affected equations (nonzero resid or grad):',len(rows))
nfail_now=sum(1 for i in affected if (i in F1))
print('  of which currently failing:',nfail_now)

# 6. build dense mod-p matrix (rows x len(cols)) and solve for consistency
nr=len(rows); nc=len(cols)
print(f'system {nr} x {nc}')
# Use flint fmpz_mod_mat. Augment with rhs, row-reduce, check consistency.
ctx=flint.fmpz_mod_ctx(p)
M=flint.fmpz_mod_mat(nr,nc+1,ctx)
for r,(g,b) in enumerate(zip(rows,rhs)):
    for c,coef in g.items():
        M[r,c]=coef
    M[r,nc]=b
R=M.rref()[0]
# consistency: any row [0...0 | nonzero]?
incons=0; rank=0
for r in range(nr):
    nzc=[c for c in range(nc) if int(R[r,c])!=0]
    if nzc: rank+=1
    elif int(R[r,nc])!=0: incons+=1
print('rank',rank,' inconsistent rows',incons)
if incons==0:
    print('TANGENT SYSTEM CONSISTENT -> mod-p Newton step exists')
else:
    print('TANGENT SYSTEM INCONSISTENT at this state (first-order)')
