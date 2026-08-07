"""URGENT: greedy-ISD at a window level, collecting EVERY code support of weight <= 6,
then testing each one: (i) mod-p consistency of the retained rows, (ii) exact HNF over Z.
Any integral one is applied and checked with solve_lab/checker.py."""
import sys, json, time, random, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
path=sys.argv[1]; LEV=int(sys.argv[2]); WMAX=int(sys.argv[3]); TL=float(sys.argv[4])
v=L.load(path); fe=L.failing_eqs(L.all_atom_values(v)); s0=L.NEQ-len(fe)
A=set(a for e in fe for a in L.eq_atoms[e][2])
for _ in range(LEV):
    R=set()
    for a in A: R|=set(L.atom2eq[a])
    A=set(a for e in R for a in L.eq_atoms[e][2])
K,Rr,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
assert len(good)==len(rows)
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
Nq=[[x%Q for x in r] for r in N]; Np=[[x%P for x in r] for r in N]; Bp=[x%P for x in B]
print('lev%d n=%d nk=%d (state score %d)'%(LEV,n,nk,s0),flush=True)
def consist_modp(idx):
    M=[[Np[i][j] for j in range(nk)]+[Bp[i]] for i in idx]; r=0
    for c in range(nk):
        pr=None
        for i in range(r,len(M)):
            if M[i][c]: pr=i;break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        inv=pow(M[r][c],-1,P); M[r]=[x*inv%P for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                f=M[i][c]; M[i]=[(a-f*b)%P for a,b in zip(M[i],M[r])]
        r+=1
    return all(M[i][nk]==0 for i in range(r,len(M)))
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
# --- direct: single-knob columns give codewords for free ---
colsup=[(sum(1 for i in range(n) if N[i][j]),j) for j in range(nk)]
colsup.sort()
print('smallest single-knob column supports: %s'%[(K[j],s) for s,j in colsup[:8]],flush=True)
cands=set()
for s,j in colsup:
    if s<=WMAX: cands.add(frozenset(i for i in range(n) if N[i][j]))
def greedy(order):
    T=[]; I=[]; piv=[]
    for i in order:
        w=Nq[i][:]
        for k2,c in enumerate(piv):
            if w[c]:
                f=w[c]; w=[(a-f*b)%Q for a,b in zip(w,T[k2])]
        nzc=[c for c in range(nk) if w[c]]
        if not nzc: continue
        c0=nzc[0]; inv=pow(w[c0],Q-2,Q); w=[x*inv%Q for x in w]
        for k2 in range(len(T)):
            if T[k2][c0]:
                f=T[k2][c0]; T[k2]=[(a-f*b)%Q for a,b in zip(T[k2],w)]
        T.append(w); piv.append(c0); I.append(i)
        if len(I)==nk: break
    return I if len(I)==nk else None
random.seed(97); t0=time.time(); tr=0; minw=10**9
order0=list(range(n))
while time.time()-t0<TL:
    tr+=1; random.shuffle(order0)
    I=greedy(order0)
    if I is None: continue
    Aug=[Nq[i][:]+[1 if k2==t else 0 for t in range(nk)] for k2,i in enumerate(I)]
    rr=0; pv=[]
    for c in range(nk):
        pr=None
        for i2 in range(rr,nk):
            if Aug[i2][c]: pr=i2;break
        if pr is None: continue
        Aug[rr],Aug[pr]=Aug[pr],Aug[rr]
        inv=pow(Aug[rr][c],Q-2,Q); Aug[rr]=[x*inv%Q for x in Aug[rr]]
        for i2 in range(nk):
            if i2!=rr and Aug[i2][c]:
                f=Aug[i2][c]; Aug[i2]=[(a-f*b)%Q for a,b in zip(Aug[i2],Aug[rr])]
        pv.append(c); rr+=1
    if rr!=nk: continue
    Uc=[[Aug[k2][nk+t] for t in range(nk)] for k2 in range(nk)]
    for t in range(nk):
        u=[Uc[j][t] for j in range(nk)]
        cw=[sum(Nq[i2][j]*u[j] for j in range(nk))%Q for i2 in range(n)]
        s=frozenset(i2 for i2 in range(n) if cw[i2])
        if s and len(s)<minw:
            minw=len(s); print('  greedy-ISD trial %d: weight %d'%(tr,minw),flush=True)
        if 0<len(s)<=WMAX: cands.add(s)
print('greedy-ISD %d trials, lightest weight seen = %d, candidate supports <= %d : %d'%(
      tr,minw,WMAX,len(cands)),flush=True)
ALL=set(range(n)); best=s0
for D in sorted(cands,key=len):
    Z=sorted(ALL-D)
    mp=consist_modp(Z)
    x=int_solve(Z) if mp else None
    print('  |D|=%d eqs=%s  modp=%s  integral=%s'%(len(D),sorted(EQ[i] for i in D),mp,x is not None),flush=True)
    if x is not None:
        w=list(v)
        for j,u in enumerate(K): w[u]=x[j]
        s2=L.NEQ-len(L.failing_eqs(L.all_atom_values(w)))
        print('  *** SCORE %d'%s2,flush=True)
        if s2>best:
            best=s2
            out='/home/user/integer_solver/solve_lab/agentA_work/A_w6_%d.json'%s2
            json.dump({str(i):str(w[i]) for i in range(L.NVARS)},open(out,'w'))
            print('  SAVED %s'%out,flush=True)
print('best from this window = %d'%best,flush=True)
