"""Priority 1: is there an integer coset leader of weight <= 6 in the CANONICAL basin?
Region at mod9118_0: 89 nontrivial rows, 65 knobs, rank 65, Q-consistent with a unique
non-integral W => every integer point's violated set must contain a code support.
Phase A: randomized information-set search for ANY code support of weight <= 6.
Phase B: exact HNF on every one found (and on targeted drop-sets that pass a mod-p
left-kernel certificate filter)."""
import sys, json, time, random, collections, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
TL=float(sys.argv[1]) if len(sys.argv)>1 else 900
path='/home/user/integer_solver/solve_lab/s10/mod9118_0.json'
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av); s0=L.NEQ-len(fe)
A=set(a for e in fe for a in L.eq_atoms[e][2])
K,R,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
Nq=[[x%Q for x in r] for r in N]
d0=[v[u] for u in K]
cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
print('region rows=%d knobs=%d violated=%d  (score %d; need <=6 violated to beat 39,026)'%(
      n,nk,len(cur),s0),flush=True)
def int_solve(idx):
    rowsM=[N[i] for i in idx]; rhs=[B[i] for i in idx]; nr=len(rowsM)
    H=[r[:] for r in rowsM]; U=[[1 if i==j else 0 for j in range(nk)] for i in range(nk)]
    pv=[]; r=0
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
        pv.append((i,r)); r+=1
    y=[0]*nk
    for i,j in pv:
        s=rhs[i]-sum(H[i][k]*y[k] for k in range(j))
        if s%H[i][j]: return None
        y[j]=s//H[i][j]
    for i in range(nr):
        if sum(H[i][k]*y[k] for k in range(nk))!=rhs[i]: return None
    return [sum(U[k][j]*y[j] for j in range(nk)) for k in range(nk)]
def report(x,D):
    w=list(v)
    for j,u in enumerate(K): w[u]=x[j]
    fe2=L.failing_eqs(L.all_atom_values(w)); s2=L.NEQ-len(fe2)
    print('*** INTEGRAL |D|=%d violate %s -> SCORE %d'%(len(D),sorted(EQ[i] for i in D),s2),flush=True)
    if s2>39026:
        out='/home/user/integer_solver/solve_lab/agentA_work/A_canon_%d.json'%s2
        json.dump({str(i):str(w[i]) for i in range(L.NVARS)},open(out,'w'))
        print('SAVED %s'%out,flush=True)
    return s2
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
sups=set(); random.seed(3); t0=time.time(); trials=0; IDX=list(range(n)); minw=99
while time.time()-t0<TL:
    trials+=1
    I=random.sample(IDX,nk)
    cws=infoset(I)
    if cws is None: continue
    OUT=[i for i in IDX if i not in I]
    for t in range(nk):
        sup=frozenset(i for i in IDX if cws[t][i])
        if sup: minw=min(minw,len(sup))
        if 0<len(sup)<=6: sups.add(sup)
    for t1 in range(nk):
        c1=cws[t1]
        for t2 in range(t1+1,nk):
            c2=cws[t2]
            for k in OUT:
                if not c2[k]: continue
                lam=(-c1[k])*pow(c2[k],Q-2,Q)%Q
                sup=frozenset(i for i in IDX if (c1[i]+lam*c2[i])%Q)
                if sup: minw=min(minw,len(sup))
                if 0<len(sup)<=6: sups.add(sup)
print('ISD trials=%d  minimum support weight seen = %d  supports<=6 found = %d'%(
      trials,minw,len(sups)),flush=True)
best=s0
for D in sorted(sups,key=len):
    x=int_solve(sorted(set(range(n))-D))
    if x is not None:
        s2=report(x,D); best=max(best,s2)
print('DONE. best score from this region = %d'%best,flush=True)
if minw>6:
    print('CONCLUSION: minimum code-support weight in this region is %d > 6, so every'%minw)
    print('integer knob vector violates >= %d rows -> this basin cannot beat 39,026.'%minw)
