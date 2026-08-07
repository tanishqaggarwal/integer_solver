"""Priority 1, decisive phase.  The mod-p necessary condition for the canonical basin is
exactly a minimum-weight syndrome-decoding problem over F_p:

    minimise |supp(x)|  subject to  Wb * x = g,      Wb : w x n over F_p.

|supp(x)| = the number of equations that must be violated.  Prange information-set
decoding: a random size-w column set I that contains the true support recovers it, with
probability ~ prod (w-k)/(n-k) per trial.  Every solution found is then HNF-tested over Z."""
import sys, json, time, random, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P
path=sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/mod9118_0.json'
TRIALS=int(sys.argv[2]) if len(sys.argv)>2 else 400
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av); s0=L.NEQ-len(fe)
A=set(a for e in fe for a in L.eq_atoms[e][2])
K,R,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
Nz=[[lin.get(j,0)%P for j in range(nk)] for e,c,lin in NZ]
Bz=[(-c)%P for e,c,lin in NZ]
Nint=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
Bint=[-c for e,c,lin in NZ]
# left kernel of N mod p
M=[[Nz[i][j] for j in range(nk)]+[1 if k==i else 0 for k in range(n)] for i in range(n)]
r=0
for c in range(nk):
    pr=None
    for i in range(r,n):
        if M[i][c]: pr=i;break
    if pr is None: continue
    M[r],M[pr]=M[pr],M[r]
    inv=pow(M[r][c],-1,P); M[r]=[x*inv%P for x in M[r]]
    for i in range(n):
        if i!=r and M[i][c]:
            f=M[i][c]; M[i]=[(a-f*b)%P for a,b in zip(M[i],M[r])]
    r+=1
Wb=[M[i][nk:] for i in range(r,n)]; w=len(Wb)
g=[sum(Wb[j][i]*Bz[i] for i in range(n))%P for j in range(w)]
cols=[[Wb[j][i] for j in range(w)] for i in range(n)]
pr_prob=1.0
for k in range(6): pr_prob*= (w-k)/(n-k)
print('%s: n=%d rows, nk=%d knobs, rank(N mod p)=%d, syndrome space w=%d'%(
      path.split('/')[-1],n,nk,r,w),flush=True)
print('Prange detection probability for a weight-6 solution: %.4f per trial'%pr_prob,flush=True)
def solve_on(I):
    """solve sum_{i in I} x_i col_i = g exactly on F_p; returns x (dict) or None"""
    A_=[[cols[i][t] for i in I]+[g[t]] for t in range(w)]
    m_=len(I); piv=[]; rr=0
    for c in range(m_):
        prw=None
        for i in range(rr,w):
            if A_[i][c]: prw=i;break
        if prw is None: continue
        A_[rr],A_[prw]=A_[prw],A_[rr]
        inv=pow(A_[rr][c],-1,P); A_[rr]=[x*inv%P for x in A_[rr]]
        for i in range(w):
            if i!=rr and A_[i][c]:
                f=A_[i][c]; A_[i]=[(a-f*b)%P for a,b in zip(A_[i],A_[rr])]
        piv.append(c); rr+=1
    for i in range(rr,w):
        if A_[i][m_]: return None
    x={}
    for i,c in enumerate(piv):
        if A_[i][m_]: x[I[c]]=A_[i][m_]
    return x
def int_solve(idx):
    rowsM=[Nint[i] for i in idx]; rhs=[Bint[i] for i in idx]; nr=len(rowsM)
    H=[q[:] for q in rowsM]; U=[[1 if a==b else 0 for b in range(nk)] for a in range(nk)]
    pv=[]; rr=0
    for i in range(nr):
        if rr>=nk: break
        while True:
            nzc=[j for j in range(rr,nk) if H[i][j]]
            if len(nzc)<=1: break
            nzc.sort(key=lambda j: abs(H[i][j])); j0=nzc[0]
            for j in nzc[1:]:
                q=H[i][j]//H[i][j0]
                if q:
                    for k2 in range(nr): H[k2][j]-=q*H[k2][j0]
                    for k2 in range(nk): U[k2][j]-=q*U[k2][j0]
        nzc=[j for j in range(rr,nk) if H[i][j]]
        if not nzc: continue
        j0=nzc[0]
        if j0!=rr:
            for k2 in range(nr): H[k2][rr],H[k2][j0]=H[k2][j0],H[k2][rr]
            for k2 in range(nk): U[k2][rr],U[k2][j0]=U[k2][j0],U[k2][rr]
        pv.append((i,rr)); rr+=1
    y=[0]*nk
    for i,j in pv:
        s=rhs[i]-sum(H[i][k2]*y[k2] for k2 in range(j))
        if s%H[i][j]: return None
        y[j]=s//H[i][j]
    for i in range(nr):
        if sum(H[i][k2]*y[k2] for k2 in range(nk))!=rhs[i]: return None
    return [sum(U[k2][j]*y[j] for j in range(nk)) for k2 in range(nk)]
random.seed(17); best=99; bestS=None; t0=time.time(); ok=0
IDX=list(range(n)); seen=set()
for tr in range(TRIALS):
    I=sorted(random.sample(IDX,w))
    x=solve_on(I)
    if x is None: continue
    ok+=1
    S=frozenset(x)
    if len(S)<best:
        best=len(S); bestS=S
        print('  trial %d: weight %d support %s [%.0fs]'%(tr,best,sorted(EQ[i] for i in S),time.time()-t0),flush=True)
    if 0<len(S)<=6 and S not in seen:
        seen.add(S)
        xi=int_solve(sorted(set(range(n))-S))
        if xi is not None:
            wv=list(v)
            for j,u in enumerate(K): wv[u]=xi[j]
            s2=L.NEQ-len(L.failing_eqs(L.all_atom_values(wv)))
            print('*** INTEGRAL |D|=%d -> SCORE %d'%(len(S),s2),flush=True)
            if s2>39026:
                out='/home/user/integer_solver/solve_lab/agentA_work/A_canon_%d.json'%s2
                json.dump({str(i):str(wv[i]) for i in range(L.NVARS)},open(out,'w'))
                print('SAVED %s'%out,flush=True); sys.exit(0)
print('trials=%d solvable=%d  MINIMUM WEIGHT OBSERVED = %d  support=%s  [%.0fs]'%(
      TRIALS,ok,best,sorted(EQ[i] for i in bestS) if bestS else None,time.time()-t0),flush=True)
miss=(1-pr_prob)**ok
print('P(a weight-6 solution exists and was missed) <= %.2g'%miss,flush=True)
# NOTE: `best` is the smallest weight SEEN and is an UPPER bound on the true minimum.
# The sound inference is: no weight-<=6 solution was found in `ok` successful trials, and
# P(one exists and was missed) <= (1-pr_prob)^ok, printed above.
