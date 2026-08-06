#!/usr/bin/env python3
"""Extract the exact conserved residue: combinations of the 16 leaf-ripple equations
whose gradient w.r.t. EVERY free input is zero, yet whose residual is nonzero.
That conserved value is the true obstruction. Then test what (if anything) moves it."""
import heal_harness as H
from jac_lib import D
import flint
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[17325]=0; H.val[9413]=0; H.forward()
H.val[7068]=H.val[2099]; H.val[4432]=H.val[19964]; H.forward()
leaf=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]

# gradients of the 16 w.r.t. ALL free inputs (seed all)
freelist=sorted(H.freeinp); colidx={j:k for k,j in enumerate(freelist)}
val=H.val
vd=[None]*H.NVARS
for j in H.freeinp: vd[j]=D(val[j],{colidx[j]:1})
ns={'v':vd,'__builtins__':{}}
for k,t in enumerate(H.order):
    r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
R=[]; G=[]
for i in leaf:
    rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
    if isinstance(rr,D): R.append(rr.v); G.append(rr.g)
    else: R.append(rr%p); G.append({})
# support = union of gradient columns
supp=sorted(set().union(*[set(g) for g in G]))
print('16 leaf-ripple: gradient support (distinct free cols):',len(supp))
sidx={c:k for k,c in enumerate(supp)}
ns_=len(supp)
# Build 16 x ns_ matrix J and find LEFT null space (combos c: c^T J = 0)
ctx=flint.fmpz_mod_ctx(p)
# left null space of J = null space of J^T. Build Jt (ns_ x 16), rref, get null.
# Easier: augment [J | R] (16 x ns_+1), row reduce; any pivot-free structure in R reveals conserved.
# We want c (1x16) with c J=0 and cR!=0. Compute rank(J) and rank([J|R]).
J=flint.fmpz_mod_mat(16,ns_,ctx)
for r,g in enumerate(G):
    for c,co in g.items(): J[r,sidx[c]]=co%p
JR=flint.fmpz_mod_mat(16,ns_+1,ctx)
for r,g in enumerate(G):
    for c,co in g.items(): JR[r,sidx[c]]=co%p
    JR[r,ns_]=R[r]%p
rJ=J.rank(); rJR=JR.rank()
print(f'rank(J)={rJ if False else rJ}' if False else f'rank(J)={rJ}'.replace('rJ',str(rJ)))
print(f'rank([J|R])={rJR}')
print('conserved-residue dimension (rank[J|R]-rank[J]):', rJR-rJ, '(1 => single conserved obstruction)')
# dimension of left-null of J that hits R:
print('left-null(J) dim =', 16-rJ, ';  of which annihilate R too =', 16-rJR)
# So #conserved that are OBSTRUCTIONS (nonzero on R) = (16-rJ) - (16-rJR) = rJR-rJ
