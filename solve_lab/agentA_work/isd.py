"""Exact max-satisfy on the residual region by information-set decoding over Z."""
import sys, json, random, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from agrow import model
import amk_model as MK
P=env.P
EXTRA=[37887,41906]
A,K,R,rows,QUAD=model(EXTRA)
aff=[(e,c,lin) for e,c,lin,hq in rows if not hq]
nk=len(K); n=len(aff)
EQ=[e for e,_,_ in aff]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in aff]
B=[-c for e,c,lin in aff]
d0=[MK.v0[u] for u in K]
cur=[sum(N[i][j]*d0[j] for j in range(nk))-B[i] for i in range(n)]
base=sum(1 for x in cur if x==0)
print('rows=%d knobs=%d baseline satisfied=%d violated=%s'%(n,nk,base,[EQ[i] for i in range(n) if cur[i]]),flush=True)

def int_solve(rowsM, rhs, nk):
    nr=len(rowsM)
    H=[r[:] for r in rowsM]
    U=[[1 if i==j else 0 for j in range(nk)] for i in range(nk)]
    piv=[]; r=0
    for i in range(nr):
        if r>=nk: break
        while True:
            nz=[j for j in range(r,nk) if H[i][j]]
            if len(nz)<=1: break
            nz.sort(key=lambda j: abs(H[i][j]))
            j0=nz[0]
            for j in nz[1:]:
                q=H[i][j]//H[i][j0]
                if q:
                    for k in range(nr): H[k][j]-=q*H[k][j0]
                    for k in range(nk): U[k][j]-=q*U[k][j0]
        nz=[j for j in range(r,nk) if H[i][j]]
        if not nz: continue
        j0=nz[0]
        if j0!=r:
            for k in range(nr): H[k][r],H[k][j0]=H[k][j0],H[k][r]
            for k in range(nk): U[k][r],U[k][j0]=U[k][j0],U[k][r]
        piv.append((i,r)); r+=1
    y=[0]*nk
    for i,j in piv:
        s=rhs[i]-sum(H[i][k]*y[k] for k in range(j))
        if s % H[i][j]: return None
        y[j]=s//H[i][j]
    for i in range(nr):
        if sum(H[i][k]*y[k] for k in range(nk))!=rhs[i]: return None
    return [sum(U[k][j]*y[j] for j in range(nk)) for k in range(nk)]

def viol(x):
    return [EQ[i] for i in range(n) if sum(N[i][j]*x[j] for j in range(nk))!=B[i]]

best=base; bx=None
t0=time.time(); tries=0; solved=0
random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 1)
IDX=list(range(n))
LIM=float(sys.argv[1]) if len(sys.argv)>1 else 600
while time.time()-t0 < LIM:
    tries+=1
    S=random.sample(IDX, nk)
    x=int_solve([N[i] for i in S],[B[i] for i in S],nk)
    if x is None: continue
    solved+=1
    vv=viol(x); s=n-len(vv)
    if s>best:
        best=s; bx=x
        print('  new best %d satisfied, violated=%s (try %d)'%(s,vv,tries),flush=True)
        json.dump({'K':K,'x':[str(t) for t in x],'sat':s},
                  open('/home/user/integer_solver/solve_lab/agentA_work/isd_best.json','w'))
print('tries=%d int-solvable=%d best=%d/%d (baseline %d)'%(tries,solved,best,n,base),flush=True)
