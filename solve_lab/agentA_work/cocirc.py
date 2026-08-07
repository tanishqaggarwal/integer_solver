"""Exact max-satisfy: enumerate rank-21 flats (hyperplanes) of the 40 region rows,
i.e. cocircuits D = complement.  For each small D test whether the affine solution
space W + ker(N_{V\D}) contains an INTEGER point (exact HNF)."""
import sys, json, random, time, collections, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from agrow import model
import amk_model as MK
P=env.P
A,K,R,rows,QUAD=model([37887,41906])
aff=[(e,c,lin) for e,c,lin,hq in rows if not hq]
nk=len(K); n=len(aff)
EQ=[e for e,_,_ in aff]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in aff]
B=[-c for e,c,lin in aff]
Q=(1<<61)-1
Nq=[[x%Q for x in r] for r in N]
d0=[MK.v0[u] for u in K]
cur=set(EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i])
print('rows=%d knobs=%d currently violated %s'%(n,nk,sorted(cur)),flush=True)

def rank_and_basis(idx):
    M=[Nq[i][:] for i in idx]; r=0; piv=[]
    for c in range(nk):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],Q-2,Q); M[r]=[x*inv%Q for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%Q for a,b in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==len(M): break
    return r,M[:r],piv

def in_span(vec,Mr,piv):
    w=vec[:]
    for i,c in enumerate(piv):
        if w[c]:
            f=w[c]; w=[(a-f*b)%Q for a,b in zip(w,Mr[i])]
    return all(x==0 for x in w)

def int_solve(rowsM,rhs):
    nr=len(rowsM)
    H=[r[:] for r in rowsM]
    U=[[1 if i==j else 0 for j in range(nk)] for i in range(nk)]
    piv=[]; r=0
    for i in range(nr):
        if r>=nk: break
        while True:
            nz=[j for j in range(r,nk) if H[i][j]]
            if len(nz)<=1: break
            nz.sort(key=lambda j: abs(H[i][j])); j0=nz[0]
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
        if s%H[i][j]: return None
        y[j]=s//H[i][j]
    for i in range(nr):
        if sum(H[i][k]*y[k] for k in range(nk))!=rhs[i]: return None
    return [sum(U[k][j]*y[j] for j in range(nk)) for k in range(nk)]

seen=set(); found=[]
random.seed(7)
t0=time.time()
LIM=float(sys.argv[1]) if len(sys.argv)>1 else 300
IDX=list(range(n)); it=0
while time.time()-t0<LIM:
    it+=1
    S=random.sample(IDX,21)
    r,Mr,piv=rank_and_basis(S)
    if r!=21: continue
    flat=[i for i in range(n) if in_span(Nq[i],Mr,piv)]
    D=tuple(sorted(set(range(n))-set(flat)))
    if D in seen: continue
    seen.add(D)
    found.append(D)
print('sampled %d, distinct cocircuits %d, sizes %s'%(it,len(found),
      sorted(collections.Counter(len(d) for d in found).items())),flush=True)
small=sorted([d for d in found if len(d)<=8],key=len)
print('cocircuits of size <=8: %d'%len(small))
best=n-len(cur)
for D in small:
    Z=[i for i in range(n) if i not in D]
    x=int_solve([N[i] for i in Z],[B[i] for i in Z])
    ok='INTEGRAL' if x is not None else 'no'
    print('   |D|=%d  violate %s  -> %s'%(len(D),[EQ[i] for i in D],ok),flush=True)
    if x is not None and n-len(D)>best:
        best=n-len(D)
        json.dump({'K':K,'x':[str(t) for t in x],'violate':[EQ[i] for i in D]},
                  open('/home/user/integer_solver/solve_lab/agentA_work/cocirc_best.json','w'))
        print('   *** IMPROVEMENT: satisfies %d of %d ***'%(best,n),flush=True)
print('best satisfied %d of %d (baseline %d)'%(best,n,n-len(cur)))
json.dump([list(d) for d in found],open('/home/user/integer_solver/solve_lab/agentA_work/cocircuits.json','w'))
