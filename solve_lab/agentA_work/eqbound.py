"""EQUATION-LEVEL bound, off manifold.
Window = equation-closure level L around a state: atoms of the failing equations, then
repeatedly A := atoms(eqs(A)).  Every atom in the window is FREE to be nonzero and to
cancel against the others; rows are exact affine forms for the EQUATION VALUES.
Lemma: rank(N) = #knobs and Q-consistency with a unique non-integral solution W force
every integer point's violated set D to contain a code support.
Two independent necessary conditions on D, each with an honest miss probability:
  (1) code support   -- information-set decoding on the equation-level code;
  (2) mod-p consistency -- V\\D consistent mod p iff g in span{col_i(Wb) : i in D},
      a minimum-weight syndrome-decoding problem over F_p (exhaustive for small |D|,
      Prange above that).
Anything surviving both is HNF-tested over Z and, if integral, applied and checked."""
import sys, json, time, random, itertools, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
path=sys.argv[1]; LEV=int(sys.argv[2]); WMAX=int(sys.argv[3]) if len(sys.argv)>3 else 6
TRIALS=int(sys.argv[4]) if len(sys.argv)>4 else 300
EXH=int(sys.argv[5]) if len(sys.argv)>5 else 3
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av); s0=L.NEQ-len(fe)
A=set(a for e in fe for a in L.eq_atoms[e][2])
for lev in range(LEV):
    R=sorted(set(e for a in A for e in L.atom2eq[a]))
    A=set(a for e in R for a in L.eq_atoms[e][2])
K,Rr,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
assert len(good)==len(rows), 'window not exactly affine'
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
d0=[v[u] for u in K]
cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
print('%s lev%d: atoms=%d eqs=%d nontrivial rows=%d knobs=%d violated=%d (score %d)'%(
      path.split('/')[-1],LEV,len(A),len(Rr),n,nk,len(cur),s0),flush=True)
Nq=[[x%Q for x in r] for r in N]; Np=[[x%P for x in r] for r in N]; Bp=[x%P for x in B]

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

def check(D):
    x=int_solve(sorted(set(range(n))-set(D)))
    if x is None: return None
    w=list(v)
    for j,u in enumerate(K): w[u]=x[j]
    s2=L.NEQ-len(L.failing_eqs(L.all_atom_values(w)))
    print('*** INTEGRAL |D|=%d violate %s -> SCORE %d'%(len(D),sorted(EQ[i] for i in D),s2),flush=True)
    if s2>39026:
        out='/home/user/integer_solver/solve_lab/agentA_work/A_eqwin_%d.json'%s2
        json.dump({str(i):str(w[i]) for i in range(L.NVARS)},open(out,'w'))
        print('SAVED %s'%out,flush=True)
    return s2

# ---------- condition 2: mod-p ----------
M=[[Np[i][j] for j in range(nk)]+[1 if k==i else 0 for k in range(n)] for i in range(n)]
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
g=[sum(Wb[j][i]*Bp[i] for i in range(n))%P for j in range(w)]
cols=[[Wb[j][i] for j in range(w)] for i in range(n)]
zerop=sum(1 for j in range(nk) if all(Np[i][j]==0 for i in range(n)))
print('mod p: rank=%d of %d knobs (%d knob columns vanish), syndrome dim w=%d, g nonzero=%s'%(
      r,nk,zerop,w,not all(x==0 for x in g)),flush=True)
def in_span(S):
    T=[cols[i][:] for i in S]; rr=0; piv=[]
    for c in range(w):
        pr=None
        for i in range(rr,len(T)):
            if T[i][c]: pr=i;break
        if pr is None: continue
        T[rr],T[pr]=T[pr],T[rr]
        inv=pow(T[rr][c],-1,P); T[rr]=[x*inv%P for x in T[rr]]
        for i in range(len(T)):
            if i!=rr and T[i][c]:
                f=T[i][c]; T[i]=[(a-f*b)%P for a,b in zip(T[i],T[rr])]
        piv.append(c); rr+=1
        if rr==len(S): break
    y=g[:]
    for i,c in enumerate(piv):
        if y[c]:
            f=y[c]; y=[(a-f*b)%P for a,b in zip(y,T[i])]
    return all(x==0 for x in y)
if all(x==0 for x in g):
    print('  (g == 0: the mod-p filter is vacuous in this window)',flush=True)
else:
    t0=time.time(); hit=None
    for k in range(1,EXH+1):
        cnt=0
        for S in itertools.combinations(range(n),k):
            cnt+=1
            if in_span(S): hit=S; break
        print('  mod-p exhaustive |D|=%d: %d subsets, consistent: %s [%.0fs]'%(
            k,cnt,sorted(EQ[i] for i in hit) if hit else 'NONE',time.time()-t0),flush=True)
        if hit: break
    if not hit:
        print('  RIGOROUS: no |D| <= %d is mod-p consistent -> every integer point violates >= %d'%(EXH,EXH+1),flush=True)
    # Prange on the syndrome problem
    pr_prob=1.0
    for k in range(WMAX): pr_prob*=(w-k)/(n-k)
    random.seed(23); okc=0; best=10**9; t0=time.time()
    for tr in range(TRIALS):
        I=sorted(random.sample(range(n),w))
        Aa=[[cols[i][t] for i in I]+[g[t]] for t in range(w)]
        m_=len(I); piv=[]; rr=0
        for c in range(m_):
            prw=None
            for i in range(rr,w):
                if Aa[i][c]: prw=i;break
            if prw is None: continue
            Aa[rr],Aa[prw]=Aa[prw],Aa[rr]
            inv=pow(Aa[rr][c],-1,P); Aa[rr]=[x*inv%P for x in Aa[rr]]
            for i in range(w):
                if i!=rr and Aa[i][c]:
                    f=Aa[i][c]; Aa[i]=[(a-f*b)%P for a,b in zip(Aa[i],Aa[rr])]
            piv.append(c); rr+=1
        if any(Aa[i][m_] for i in range(rr,w)): continue
        okc+=1
        S=[I[c] for i,c in enumerate(piv) if Aa[i][m_]]
        if len(S)<best:
            best=len(S)
            print('  Prange trial %d: weight %d [%.0fs]'%(tr,best,time.time()-t0),flush=True)
            if best<=WMAX: check(S)
    print('  Prange: %d trials, %d solvable, lightest weight SEEN = %d (an upper bound)'%(TRIALS,okc,best),flush=True)
    print('  P(a weight-<=%d mod-p-consistent D exists and was missed) <= %.3g'%(WMAX,(1-pr_prob)**okc),flush=True)
