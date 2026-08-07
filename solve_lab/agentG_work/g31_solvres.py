"""Solve the small residual system in its own variables over F_p."""
import os, sys, pickle, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gsym2 as G
from gsym2 import L, ad, P
d=pickle.load(open(sys.argv[1],'rb')); res=d['res']; NB=d['NB']
vars_=sorted({NB[k] for a,g in res if not isinstance(g,int) for m in g for k,_ in m})
ix={u:i for i,u in enumerate(vars_)}
print('vars',vars_)
polys=[(a,g) for a,g in res if not isinstance(g,int)]
consts=[(a,g%P) for a,g in res if isinstance(g,int) and g%P]
if consts: print('UNCONDITIONAL NONZERO CONSTANTS:',consts)
lin=[(a,g) for a,g in polys if G.deg(g)==1]
print('linear residual eqs:',[a for a,_ in lin])
# solve linear subsystem
import itertools
nv=len(vars_)
rowsM=[]
for a,g in lin:
    row=[0]*(nv+1)
    for m,c in g.items():
        if not m: row[nv]=(-c)%P
        else: row[ix[NB[m[0][0]]]]=c%P
    rowsM.append(row)
def rref(M,nc):
    M=[r[:] for r in M]; piv=[]; r=0
    for c in range(nc):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]%P: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        iv=pow(M[r][c],-1,P); M[r]=[x*iv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%P:
                f=M[i][c]; M[i]=[(x-f*y)%P for x,y in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==len(M): break
    return M,piv,r
if rowsM:
    M,piv,rk=rref(rowsM,nv)
    inc=[i for i in range(len(M)) if all(x%P==0 for x in M[i][:nv]) and M[i][nv]%P]
    print('linear rank',rk,'inconsistent',len(inc))
    if inc: print('LINEAR RESIDUAL INCONSISTENT'); sys.exit()
    if rk==nv:
        sol=[0]*nv
        for r_,c in enumerate(piv): sol[c]=M[r_][nv]%P
        print('unique candidate:',{vars_[i]:sol[i] for i in range(nv)})
        ok=True
        for a,g in polys:
            val=0
            for m,c in g.items():
                t=c
                for k,e in m: t=t*pow(sol[ix[NB[k]]],e,P)%P
                val=(val+t)%P
            print('   a%-6d -> %s'%(a,'ZERO' if val==0 else val))
            if val: ok=False
        print('ALL RESIDUAL CHECKS ZERO' if ok else 'NOT a solution')
        pickle.dump({'sol':{vars_[i]:sol[i] for i in range(nv)},'flip':d['flip']},open(sys.argv[1].replace('res_','sol_'),'wb'))
    else:
        print('linear rank %d < %d : need to solve the nonlinear part'%(rk,nv))
