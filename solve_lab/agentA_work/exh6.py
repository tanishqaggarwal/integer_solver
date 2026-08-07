"""Exhaustive proof-or-refutation: is there an integer knob vector violating <= 6 of the
30 nontrivial region rows?  Enumerate every admissible violation set D (|D|<=6, i.e. D
must contain a code support), filter by mod-p consistency, then exact HNF."""
import sys, json, itertools, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from agrow import model
import amk_model as MK
P=env.P
A,K,R,rows,QUAD=model([37887,41906])
aff=[(e,c,lin) for e,c,lin,hq in rows if not hq]
nk=len(K)
NZ=[(e,c,lin) for e,c,lin in aff if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
E2I={e:i for i,e in enumerate(EQ)}
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
Np=[[x%P for x in r] for r in N]; Bp=[x%P for x in B]
sups=[frozenset(E2I[e] for e in s) for s in json.load(open('/home/user/integer_solver/solve_lab/agentA_work/supports.json'))]
print('rows=%d knobs=%d code supports loaded=%d (sizes %s)'%(n,nk,len(sups),
      sorted(collections.Counter(len(s) for s in sups).items())),flush=True)

def consist_modp(idx):
    M=[[Np[i][j] for j in range(nk)]+[Bp[i]] for i in idx]
    r=0
    for c in range(nk):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],P-2,P); M[r]=[x*inv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%P for a,b in zip(M[i],M[r])]
        r+=1
    return all(M[i][nk]==0 for i in range(r,len(M)))

def int_solve(idx):
    rowsM=[N[i] for i in idx]; rhs=[B[i] for i in idx]; nr=len(rowsM)
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

ALL=set(range(n))
cands=set()
for s in sups:
    if len(s)<=6: cands.add(s)
    if len(s)==5:
        for i in ALL-s: cands.add(s|{i})
print('candidate violation sets of size <=6: %d'%len(cands),flush=True)
t0=time.time(); passed=0; hits=[]
for k,D in enumerate(sorted(cands,key=len)):
    Z=sorted(ALL-D)
    if not consist_modp(Z): continue
    passed+=1
    x=int_solve(Z)
    if x is not None:
        hits.append((len(D),sorted(EQ[i] for i in D),x))
        print('*** INTEGRAL with |D|=%d  violate %s'%(len(D),sorted(EQ[i] for i in D)),flush=True)
        json.dump({'K':K,'x':[str(t) for t in x],'violate':sorted(EQ[i] for i in D)},
                  open('/home/user/integer_solver/solve_lab/agentA_work/exh6_hit.json','w'))
        break
    if k%2000==0: print('  %d/%d checked, %d passed mod-p [%.0fs]'%(k,len(cands),passed,time.time()-t0),flush=True)
print('DONE: %d candidates, %d passed mod-p filter, %d integral'%(len(cands),passed,len(hits)),flush=True)
