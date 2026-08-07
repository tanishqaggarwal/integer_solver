"""DEFINITIVE region computation: 31 affine rows (30 nontrivial equation rows + the
linear form Q with a37887 = Q^2, i.e. eq 8680), 22 knobs.
(1) rank / Q-consistency ; (2) EXHAUSTIVE code supports of weight <= WMAX ;
(3) exact HNF integrality of every admissible violation set of size <= WMAX."""
import sys, json, itertools, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
from agrow import model
import amk_model as MK
P=env.P; Q=(1<<61)-1
WMAX=int(sys.argv[1]) if len(sys.argv)>1 else 6
A,K,R,rows,QUAD=model([37887,41906])
aff=[(e,c,lin) for e,c,lin,hq in rows if not hq]
nk=len(K)
qr=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/qrow.json'))
assert qr['K']==K
Qlin={i:int(x) for i,x in enumerate(qr['Q']) if int(x)}
aff.append((8680,int(qr['Qc']),Qlin))
NZ=[(e,c,lin) for e,c,lin in aff if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
d0=[MK.v0[u] for u in K]
cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
print('rows=%d knobs=%d currently violated=%d %s'%(n,nk,len(cur),cur),flush=True)
# Q-consistency / rank
def rankQ(idx,aug=False):
    M=[[F(N[i][j]) for j in range(nk)]+([F(B[i])] if aug else []) for i in idx]
    r=0; nc=nk+(1 if aug else 0)
    for c in range(nk):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]!=0: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        pv=M[r][c]; M[r]=[x/pv for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                f=M[i][c]; M[i]=[x-f*y for x,y in zip(M[i],M[r])]
        r+=1
    if aug: return r, [i for i in range(r,len(M)) if M[i][nk]!=0]
    return r
r0=rankQ(range(n)); r1,inc=rankQ(range(n),aug=True)
print('rank(N)=%d  inconsistent rows over Q: %d'%(r0,len(inc)),flush=True)
Nq=[[x%Q for x in r] for r in N]
# information set
I=[]; M=[]
for i in range(n):
    T=[r[:] for r in M+[Nq[i][:]]]; rr=0
    for c in range(nk):
        pr=None
        for t in range(rr,len(T)):
            if T[t][c]: pr=t;break
        if pr is None: continue
        T[rr],T[pr]=T[pr],T[rr]
        inv=pow(T[rr][c],Q-2,Q); T[rr]=[x*inv%Q for x in T[rr]]
        for t in range(len(T)):
            if t!=rr and T[t][c]:
                f=T[t][c]; T[t]=[(a-f*b)%Q for a,b in zip(T[t],T[rr])]
        rr+=1
    if rr==len(M)+1: I.append(i); M.append(Nq[i][:])
    if len(I)==nk: break
OUT=[i for i in range(n) if i not in I]
print('information set %d, OUT %d'%(len(I),len(OUT)),flush=True)
aug=[M[k]+[1 if k==t else 0 for t in range(nk)] for k in range(nk)]
rr=0
for c in range(nk):
    pr=None
    for t in range(rr,nk):
        if aug[t][c]: pr=t;break
    aug[rr],aug[pr]=aug[pr],aug[rr]
    inv=pow(aug[rr][c],Q-2,Q); aug[rr]=[x*inv%Q for x in aug[rr]]
    for t in range(nk):
        if t!=rr and aug[t][c]:
            f=aug[t][c]; aug[t]=[(a-f*b)%Q for a,b in zip(aug[t],aug[rr])]
    rr+=1
Uc=[[aug[k][nk+t] for t in range(nk)] for k in range(nk)]
CW=[[sum(Nq[i][j]*Uc[j][t] for j in range(nk))%Q for i in range(n)] for t in range(nk)]
Aout=[[CW[t][k] for k in OUT] for t in range(nk)]
nout=len(OUT)
def kernel(rowsM, ncols):
    T=[r[:] for r in rowsM]; nr=len(T); piv=[]; rr=0
    for c in range(ncols):
        pr=None
        for t in range(rr,nr):
            if T[t][c]: pr=t;break
        if pr is None: continue
        T[rr],T[pr]=T[pr],T[rr]
        inv=pow(T[rr][c],Q-2,Q); T[rr]=[x*inv%Q for x in T[rr]]
        for t in range(nr):
            if t!=rr and T[t][c]:
                f=T[t][c]; T[t]=[(a-f*b)%Q for a,b in zip(T[t],T[rr])]
        piv.append(c); rr+=1
    free=[c for c in range(ncols) if c not in piv]
    out=[]
    for fc in free:
        vec=[0]*ncols; vec[fc]=1
        for k,pc in enumerate(piv): vec[pc]=(-T[k][fc])%Q
        out.append(vec)
    return out
sups=set(); t0=time.time()
for j in range(1,WMAX+1):
    need=nout-(WMAX-j)
    if need<0: need=0
    if need>nout: continue
    for J in itertools.combinations(range(nk), j):
        AJ=[[Aout[t][k] for k in range(nout)] for t in J]
        for S in itertools.combinations(range(nout), need):
            sub=[[AJ[a][s] for a in range(j)] for s in S]
            for u in kernel(sub, j):
                if all(x==0 for x in u): continue
                cw=[0]*n
                for a,t in enumerate(J):
                    if u[a]:
                        ua=u[a]
                        for i in range(n): cw[i]=(cw[i]+ua*CW[t][i])%Q
                sup=frozenset(EQ[i] for i in range(n) if cw[i])
                if 0<len(sup)<=WMAX: sups.add(sup)
    print('  j=%d supports=%d [%.0fs]'%(j,len(sups),time.time()-t0),flush=True)
print('EXHAUSTIVE supports <=%d : %d  sizes=%s'%(WMAX,len(sups),
      sorted(collections.Counter(len(s) for s in sups).items())),flush=True)
json.dump([sorted(s) for s in sups],open('/home/user/integer_solver/solve_lab/agentA_work/supports31.json','w'))
# --- integrality ---
E2I={e:i for i,e in enumerate(EQ)}
Np=[[x%P for x in r] for r in N]; Bp=[x%P for x in B]
def consist_modp(idx):
    M=[[Np[i][j] for j in range(nk)]+[Bp[i]] for i in idx]; r=0
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
ALL=set(range(n)); cands=set()
for s in sups:
    si=frozenset(E2I[e] for e in s)
    if len(si)<=WMAX: cands.add(si)
    if len(si)<WMAX:
        for extra in itertools.combinations(sorted(ALL-si), WMAX-len(si)):
            cands.add(si|set(extra))
print('admissible violation sets of size <= %d : %d'%(WMAX,len(cands)),flush=True)
t0=time.time(); passed=0
for k,D in enumerate(sorted(cands,key=len)):
    Z=sorted(ALL-D)
    if not consist_modp(Z): continue
    passed+=1
    x=int_solve(Z)
    if x is not None:
        print('*** INTEGRAL |D|=%d violate %s'%(len(D),sorted(EQ[i] for i in D)),flush=True)
        json.dump({'K':K,'x':[str(t) for t in x],'violate':sorted(EQ[i] for i in D)},
                  open('/home/user/integer_solver/solve_lab/agentA_work/full31_hit.json','w'))
        sys.exit(0)
    if k%5000==0: print('  %d/%d  passed-modp=%d [%.0fs]'%(k,len(cands),passed,time.time()-t0),flush=True)
print('DONE: %d candidates, %d passed mod-p, 0 integral -> 7 violated is OPTIMAL in this region'%(len(cands),passed),flush=True)
