"""Exhaustive-ish minimum-support search in the column space of the region matrix,
then EXACT integer solvability of the complementary row system (HNF).
Improvement over the deliverable needs a support of size <= 6 that is Z-solvable."""
import sys, json, random, time, collections, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
from agrow import model
import amk_model as MK
P=env.P
A,K,R,rows,QUAD=model([37887,41906])
aff=[(e,c,lin) for e,c,lin,hq in rows if not hq]
nk=len(K)
NZ=[(e,c,lin) for e,c,lin in aff if lin]
ZERO=[(e,c) for e,c,lin in aff if not lin]
print('nontrivial rows=%d  all-zero rows=%d (const nonzero: %s)'%(
      len(NZ),len(ZERO),[e for e,c in ZERO if c!=0]))
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
d0=[MK.v0[u] for u in K]
cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
print('currently violated: %s (%d)'%(cur,len(cur)))

def int_solve(rowsM,rhs):
    nr=len(rowsM)
    H=[r[:] for r in rowsM]; U=[[1 if i==j else 0 for j in range(nk)] for i in range(nk)]
    piv=[]; r=0
    for i in range(nr):
        if r>=nk: break
        while True:
            nzc=[j for j in range(r,nk) if H[i][j]]
            if len(nzc)<=1: break
            nzc.sort(key=lambda j: abs(H[i][j])); j0=nzc[0]
            for j in nzc[1:]:
                q=H[i][j]//H[i][j0]
                if q:
                    for k in range(nr): H[k][j]-=q*H[k][j0]
                    for k in range(nk): U[k][j]-=q*U[k][j0]
        nzc=[j for j in range(r,nk) if H[i][j]]
        if not nzc: continue
        j0=nzc[0]
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

# --- minimum-support search by random information sets over Q ---
def inv_solve(I):
    """return function mapping e_j (j in I positions) to full codeword, or None"""
    M=[[F(N[i][j]) for j in range(nk)] for i in I]
    aug=[M[k]+[F(1) if k==t else F(0) for t in range(len(I))] for k in range(len(I))]
    r=0; piv=[]
    for c in range(nk):
        pr=None
        for i in range(r,len(aug)):
            if aug[i][c]!=0: pr=i;break
        if pr is None: continue
        aug[r],aug[pr]=aug[pr],aug[r]
        pv=aug[r][c]; aug[r]=[x/pv for x in aug[r]]
        for i in range(len(aug)):
            if i!=r and aug[i][c]!=0:
                f=aug[i][c]; aug[i]=[x-f*y for x,y in zip(aug[i],aug[r])]
        piv.append(c); r+=1
    if r!=nk: return None
    # u for target vector t (values on I): u = sum over rows
    Uc=[[aug[k][nk+t] for t in range(len(I))] for k in range(nk)]   # u_j = sum_t Uc[j][t]*val_t
    return Uc

random.seed(int(sys.argv[2]) if len(sys.argv)>2 else 3)
LIM=float(sys.argv[1]) if len(sys.argv)>1 else 300
t0=time.time(); best={}
IDX=list(range(n)); trials=0
while time.time()-t0<LIM:
    trials+=1
    I=random.sample(IDX,nk)
    Uc=inv_solve(I)
    if Uc is None: continue
    OUT=[i for i in IDX if i not in I]
    cw=[]   # codewords for e_t
    for t in range(nk):
        u=[Uc[j][t] for j in range(nk)]
        c=[sum(F(N[i][j])*u[j] for j in range(nk)) for i in IDX]
        cw.append(c)
    # weight-1 in I
    for t in range(nk):
        c=cw[t]; sup=frozenset(EQ[i] for i in IDX if c[i]!=0)
        if len(sup)<=8: best[sup]=1
    # weight-2 in I, killing one outside coordinate
    for t1 in range(nk):
        for t2 in range(t1+1,nk):
            c1,c2=cw[t1],cw[t2]
            for k in OUT:
                if c2[k]==0: continue
                lam=-c1[k]/c2[k]
                c=[a+lam*b for a,b in zip(c1,c2)]
                sup=frozenset(EQ[i] for i in IDX if c[i]!=0)
                if 0<len(sup)<=8: best[sup]=2
print('trials=%d distinct low-weight supports found=%d'%(trials,len(best)),flush=True)
cands=sorted(best,key=len)
found=[]
for sup in cands:
    if len(sup)>7: continue
    D=set(sup); Z=[i for i in range(n) if EQ[i] not in D]
    x=int_solve([N[i] for i in Z],[B[i] for i in Z])
    ok = x is not None
    print('  |D|=%d %s -> %s'%(len(sup),sorted(sup),'INTEGRAL' if ok else 'no'),flush=True)
    if ok: found.append((len(sup),sorted(sup),x))
print('minimum support size found:',min(len(s) for s in cands) if cands else None)
if found:
    found.sort()
    json.dump({'K':K,'x':[str(t) for t in found[0][2]],'violate':found[0][1]},
              open('/home/user/integer_solver/solve_lab/agentA_work/minwt_best.json','w'))
    print('BEST Z-SOLVABLE: violate %d rows: %s'%(found[0][0],found[0][1]))
