import os, sys, json, pickle, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym
from gsym import *
d=pickle.load(open('sys17.pkl','rb')); SYMS=d['syms']; rows=d['rows']; n=len(SYMS)
print('symbols', SYMS)
lin=[(a,f) for a,f in rows if gsym.deg(f)==1]
non=[(a,f) for a,f in rows if gsym.deg(f)>1]
print('linear %d, nonlinear %d'%(len(lin),len(non)))
def rref(M, ncol):
    M=[r[:] for r in M]; piv=[]; r=0
    for c in range(ncol):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]%P: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        iv=pow(M[r][c],-1,P)
        M[r]=[x*iv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%P:
                f=M[i][c]; M[i]=[(x-f*y)%P for x,y in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==len(M): break
    return M,piv,r
# augmented matrix: cols 0..n-1 vars, col n = -constant
A=[]
for a,f in lin:
    row=[0]*(n+1)
    for m,c in f.items():
        if sum(m)==0: row[n]=(-c)%P
        else:
            i=[k for k,e in enumerate(m) if e][0]; row[i]=c%P
    A.append(row)
M,piv,rk=rref(A,n)
print('linear system rank', rk, 'pivots', [SYMS[c] for c in piv])
incons=[i for i in range(len(M)) if all(x%P==0 for x in M[i][:n]) and M[i][n]%P]
print('inconsistent rows:', incons)
free=[c for c in range(n) if c not in piv]
print('free params:', [SYMS[c] for c in free])
# parametrize: x_piv = M[r][n] - sum over free M[r][c]*t_c
pickle.dump({'M':M,'piv':piv,'free':free,'SYMS':SYMS,'nonlin':non}, open('linpar.pkl','wb'))
for a,f in non:
    print('  nonlinear a%-6d deg%d terms%d'%(a,gsym.deg(f),len(f)))
