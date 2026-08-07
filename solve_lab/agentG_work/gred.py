"""Reduce a given state to the exact mod-p residual: symbolic pass over the
non-boolean free inputs, linear solve, substitution into the nonlinear checks."""
import os, sys, json, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *

def rref(M,ncol):
    M=[r[:] for r in M]; piv=[]; r=0
    for c in range(ncol):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]%P: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        iv=pow(M[r][c],-1,P); M[r]=[x*iv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%P:
                fq=M[i][c]; M[i]=[(x-fq*y)%P for x,y in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==len(M): break
    return M,piv,r

def reduce_state(v, SYMS, cap=None):
    """returns dict with residual info. v must already be forward-evaluated."""
    n=len(SYMS)
    val=gsym.build(v,SYMS,cap=cap,verbose=False)
    rows=[];nzc=[]
    for a in gsym.check_atoms():
        f=gsym.evalpoly_sym(a,val,n,cap)
        if isinstance(f,int):
            if f%P: nzc.append((a,f%P))
        else: rows.append((a,f))
    lin=[(a,f) for a,f in rows if gsym.deg(f)==1]
    non=[(a,f) for a,f in rows if gsym.deg(f)>1]
    A=[]
    for a,f in lin:
        row=[0]*(n+1)
        for m,c in f.items():
            if sum(m)==0: row[n]=(-c)%P
            else:
                i=[k for k,e in enumerate(m) if e][0]; row[i]=c%P
        A.append(row)
    M,piv,rk=rref(A,n) if A else ([],[],0)
    inc=[i for i in range(len(M)) if all(x%P==0 for x in M[i][:n]) and M[i][n]%P]
    free=[c for c in range(n) if c not in piv]; k=len(free)
    sub=[None]*n
    for j,c in enumerate(free):
        e=[0]*k; e[j]=1; sub[c]={tuple(e):1}
    for r,c in enumerate(piv):
        f={}; const=M[r][n]%P
        if const: f[(0,)*k]=const
        for j,c2 in enumerate(free):
            co=(-M[r][c2])%P
            if co:
                e=[0]*k; e[j]=1; f[tuple(e)]=co
        sub[c]=f if f else 0
    def subpoly(f):
        out=0
        for m,c in f.items():
            t=c%P
            for i,e in enumerate(m):
                for _ in range(e):
                    t=gsym.pmul(t,sub[i],k,cap)
                    if t==0: break
                if t==0: break
            if t!=0: out=gsym.padd(out,t,k)
        return out
    res=[(a,subpoly(f)) for a,f in non]
    return dict(nzc=nzc,rank=rk,ninc=len(inc),nlin=len(lin),nnon=len(non),
                nfree=k,res=res,M=M,piv=piv,free=free,sub=sub,rows=rows)

def residual_signature(r):
    out=[]
    for a,g in r['res']:
        if isinstance(g,int):
            if g%P: out.append((a,g%P))
        else:
            out.append((a,'POLY%d'%len(g)))
    return out
