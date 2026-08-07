"""Necessary condition (1) at EQUATION level: information-set decoding for a code support
of weight <= WMAX in the equation-level code of a window, with an honest miss probability."""
import sys, json, time, random, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
path=sys.argv[1]; LEV=int(sys.argv[2]); WMAX=int(sys.argv[3]); TL=float(sys.argv[4])
v=L.load(path); fe=L.failing_eqs(L.all_atom_values(v)); s0=L.NEQ-len(fe)
A=set(a for e in fe for a in L.eq_atoms[e][2])
for _ in range(LEV):
    R=sorted(set(e for a in A for e in L.atom2eq[a]))
    A=set(a for e in R for a in L.eq_atoms[e][2])
K,Rr,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
assert len(good)==len(rows)
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
Nq=[[x%Q for x in r] for r in N]
pin=nk/n
pdet=sum(__import__('math').comb(WMAX,j)*pin**j*(1-pin)**(WMAX-j) for j in range(3))
print('lev%d n=%d rows nk=%d knobs ; detection prob for a weight-%d support = %.4f/trial'%(
      LEV,n,nk,WMAX,pdet),flush=True)
def int_solve(idx):
    rowsM=[N[i] for i in idx]; rhs=[B[i] for i in idx]; nr=len(rowsM)
    H=[q[:] for q in rowsM]; U=[[1 if a==b else 0 for b in range(nk)] for a in range(nk)]
    pv=[]; rr=0
    for i in range(nr):
        if rr>=nk: break
        while True:
            nzc=[j for j in range(rr,nk) if H[i][j]]
            if len(nzc)<=1: break
            nzc.sort(key=lambda j: abs(H[i][j])); j0=nzc[0]
            for j in nzc[1:]:
                q2=H[i][j]//H[i][j0]
                if q2:
                    for k2 in range(nr): H[k2][j]-=q2*H[k2][j0]
                    for k2 in range(nk): U[k2][j]-=q2*U[k2][j0]
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
def infoset(I):
    M=[[Nq[i][j] for j in range(nk)]+[1 if k==t else 0 for t in range(len(I))] for k,i in enumerate(I)]
    r=0
    for c in range(nk):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]: pr=i;break
        if pr is None: return None
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],Q-2,Q); M[r]=[x*inv%Q for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%Q for a,b in zip(M[i],M[r])]
        r+=1
    Uc=[[M[k][nk+t] for t in range(len(I))] for k in range(nk)]
    return [[sum(Nq[i][j]*Uc[j][t] for j in range(nk))%Q for i in range(n)] for t in range(nk)]
random.seed(31); t0=time.time(); tr=0; ok=0; minw=10**9; sups=set()
while time.time()-t0<TL:
    tr+=1
    I=random.sample(range(n),nk)
    cws=infoset(I)
    if cws is None: continue
    ok+=1
    OUT=[i for i in range(n) if i not in I]
    for t in range(nk):
        s=frozenset(i for i in range(n) if cws[t][i])
        if s: minw=min(minw,len(s))
        if 0<len(s)<=WMAX: sups.add(s)
    for t1 in range(nk):
        c1=cws[t1]
        for t2 in range(t1+1,nk):
            c2=cws[t2]
            for k in OUT:
                if not c2[k]: continue
                lam=(-c1[k])*pow(c2[k],Q-2,Q)%Q
                s=frozenset(i for i in range(n) if (c1[i]+lam*c2[i])%Q)
                if s: minw=min(minw,len(s))
                if 0<len(s)<=WMAX: sups.add(s)
print('trials=%d usable=%d lightest support SEEN=%d (upper bound) supports<=%d found=%d [%.0fs]'%(
      tr,ok,minw,WMAX,len(sups),time.time()-t0),flush=True)
print('P(a weight-<=%d support exists and was missed) <= %.3g'%(WMAX,(1-pdet)**ok),flush=True)
for D in sorted(sups,key=len):
    x=int_solve(sorted(set(range(n))-D))
    if x is not None:
        w=list(v)
        for j,u in enumerate(K): w[u]=x[j]
        s2=L.NEQ-len(L.failing_eqs(L.all_atom_values(w)))
        print('*** INTEGRAL |D|=%d -> SCORE %d'%(len(D),s2),flush=True)
        if s2>39026:
            out='/home/user/integer_solver/solve_lab/agentA_work/A_eqisd_%d.json'%s2
            json.dump({str(i):str(w[i]) for i in range(L.NVARS)},open(out,'w'))
            print('SAVED %s'%out,flush=True)
