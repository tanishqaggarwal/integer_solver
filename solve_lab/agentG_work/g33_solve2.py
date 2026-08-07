"""Solve a small residual system over F_p: linear elimination then univariate roots."""
import os, sys, pickle, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
import flint
d=pickle.load(open(sys.argv[1],'rb')); res=d['res']; NB=d['NB']
polys=[(a,g) for a,g in res if not isinstance(g,int)]
consts=[(a,g%P) for a,g in res if isinstance(g,int) and g%P]
print('nonzero constants:',consts)
vars_=sorted({NB[k] for a,g in polys for m in g for k,_ in m})
ix={u:i for i,u in enumerate(vars_)}; nv=len(vars_)
print('vars',vars_)
def conv(g):
    out={}
    for m,c in g.items():
        e=[0]*nv
        for k,ee in m: e[ix[NB[k]]]=ee
        out[tuple(e)]=c%P
    return out
S=[(a,conv(g)) for a,g in polys]
def deg(f): return max(sum(m) for m in f) if f else 0
def subst(f, j, expr):
    """substitute var j := expr (dict monomial->coef in the same nv-space)"""
    out={}
    for m,c in f.items():
        e=m[j]
        base=list(m); base[j]=0; base=tuple(base)
        cur={base:c}
        for _ in range(e):
            nxt={}
            for m1,c1 in cur.items():
                for m2,c2 in expr.items():
                    mm=tuple(a+b for a,b in zip(m1,m2))
                    nxt[mm]=(nxt.get(mm,0)+c1*c2)%P
            cur={k:v for k,v in nxt.items() if v}
        for k,v in cur.items():
            out[k]=(out.get(k,0)+v)%P
    return {k:v for k,v in out.items() if v}
# repeatedly: find an equation linear in some variable with constant coefficient
elim={}
changed=True
while changed:
    changed=False
    for idx,(a,f) in enumerate(S):
        if not f: continue
        for j in range(nv):
            if j in elim: continue
            # f = c*x_j + rest with c constant and x_j appearing only to power 1
            if any(m[j]>1 for m in f): continue
            cm=[m for m in f if m[j]==1]
            if not cm: continue
            if any(sum(m)!=1 for m in cm): continue   # coefficient must be constant
            c=f[cm[0]]
            rest={tuple(0 if k==j else e for k,e in enumerate(m)):v for m,v in f.items() if m[j]==0}
            inv=pow(c,-1,P)
            expr={m:(-v*inv)%P for m,v in rest.items()}
            print('  eliminate x%d using a%d'%(vars_[j],a))
            elim[j]=expr
            S=[(aa,(subst(ff,j,expr) if ff else ff)) for aa,ff in S]
            S=[(aa,ff) for aa,ff in S if aa!=a or ff]
            # drop the used equation
            S=[(aa,ff) for k,(aa,ff) in enumerate(S)]
            S[idx]=(a,{})
            changed=True
            break
        if changed: break
rem=[(a,f) for a,f in S if f]
print('remaining equations:',[(a,deg(f),len(f)) for a,f in rem])
free=[j for j in range(nv) if j not in elim]
print('free vars:',[vars_[j] for j in free])
if len(free)==1:
    j=free[0]
    ctx=flint.fmpz_mod_poly_ctx(P)
    common=None
    for a,f in rem:
        coef=[0]*(max(m[j] for m in f)+1)
        for m,c in f.items(): coef[m[j]]=(coef[m[j]]+c)%P
        pol=ctx(coef)
        common=pol if common is None else common.gcd(pol)
    print('gcd degree', common.degree())
    rts=common.roots() if common.degree()>0 else []
    print('roots:',rts)
elif len(free)==0:
    print('fully determined; residual equations left:',len(rem))
else:
    print('%d free vars left'%len(free))
pickle.dump({'S':S,'elim':elim,'vars':vars_,'free':free},open(sys.argv[1].replace('res_','sys_'),'wb'))
