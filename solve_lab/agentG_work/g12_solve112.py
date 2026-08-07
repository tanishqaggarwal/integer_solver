import os, sys, json, pickle, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
d=pickle.load(open('sys112.pkl','rb')); SYMS=d['syms']; rows=d['rows']; n=len(SYMS)
lin=[(a,f) for a,f in rows if gsym.deg(f)==1]
non=[(a,f) for a,f in rows if gsym.deg(f)>1]
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
A=[]
for a,f in lin:
    row=[0]*(n+1)
    for m,c in f.items():
        if sum(m)==0: row[n]=(-c)%P
        else:
            i=[k for k,e in enumerate(m) if e][0]; row[i]=c%P
    A.append(row)
M,piv,rk=rref(A,n)
inc=[i for i in range(len(M)) if all(x%P==0 for x in M[i][:n]) and M[i][n]%P]
print('linear: %d eqs, rank %d, inconsistent rows %d, free params %d'%(len(A),rk,len(inc),n-rk))
free=[c for c in range(n) if c not in piv]
k=len(free)
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
                t=gsym.pmul(t,sub[i],k,None)
                if t==0: break
            if t==0: break
        if t!=0: out=gsym.padd(out,t,k)
    return out
print('free params (%d): %s'%(k,[SYMS[c] for c in free]))
res=[]
for a,f in non:
    g=subpoly(f)
    res.append((a,g))
    print('a%-6d -> %s'%(a, ('CONST %d'%g) if isinstance(g,int) else 'deg %d terms %d'%(gsym.deg(g),len(g))))
pickle.dump({'res':res,'free':[SYMS[c] for c in free],'freeidx':free,'sub':sub,'SYMS':SYMS,'M':M,'piv':piv},open('red112.pkl','wb'))
