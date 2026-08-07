"""full31-style exact coset decoding on the CANONICAL basin region (mod9118_0).
Region = atoms of the failing equations; knobs strictly linear; all rows affine.
Step 1: structure (rank, consistency, non-integral coords of W).
Step 2: targeted drop-sets (equations of the residual atoms) tested by exact HNF.
Step 3: randomized ISD for any weight<=6 code support, each tested by HNF."""
import sys, json, itertools, time, collections, random
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
from fractions import Fraction as F
import env, lib as L
from regsolve2 import build, qsolve
P=env.P; Q=(1<<61)-1
path=sys.argv[1]
TL=float(sys.argv[2]) if len(sys.argv)>2 else 600
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av); s0=L.NEQ-len(fe)
nz=[a for a in range(L.NA) if av[a]]
A=set(a for e in fe for a in L.eq_atoms[e][2])
K,R,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
skipped=len(rows)-len(good)
NZR=[(e,c,lin) for e,c,lin in good if lin]
dead=[(e,c) for e,c,lin in good if not lin and c!=0]
n=len(NZR); EQ=[e for e,_,_ in NZR]; E2I={e:i for i,e in enumerate(EQ)}
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZR]
B=[-c for e,c,lin in NZR]
d0=[v[u] for u in K]
cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
print('%s score=%d failing=%d nz=%s'%(path.split('/')[-1],s0,len(fe),nz),flush=True)
print('region atoms=%d knobs=%d eqs=%d skipped=%d nontrivial=%d dead-violated=%d cur-violated=%d'%(
      len(A),nk,len(R),skipped,n,len(dead),len(cur)),flush=True)
print('violated rows:',cur,flush=True)
print('dead (unfixable) rows:',[e for e,c in dead],flush=True)
sol,free,incons,r,piv,mat,eqs=qsolve(rows,nk)
bad=[(K[j],sol[j]) for j in range(nk) if sol[j] is not None and sol[j].denominator!=1]
print('rank=%d free=%d incons=%d non-integral knobs=%s'%(r,len(free),len(incons),
      [(u,'p' if s.denominator==P else str(s.denominator)) for u,s in bad]),flush=True)

def int_solve(idx):
    rowsM=[N[i] for i in idx]; rhs=[B[i] for i in idx]; nr=len(rowsM)
    H=[rr[:] for rr in rowsM]; U=[[1 if i==j else 0 for j in range(nk)] for i in range(nk)]
    pv=[]; r2=0
    for i in range(nr):
        if r2>=nk: break
        while True:
            nzc=[j for j in range(r2,nk) if H[i][j]]
            if len(nzc)<=1: break
            nzc.sort(key=lambda j: abs(H[i][j])); j0=nzc[0]
            for j in nzc[1:]:
                q=H[i][j]//H[i][j0]
                if q:
                    for k2 in range(nr): H[k2][j]-=q*H[k2][j0]
                    for k2 in range(nk): U[k2][j]-=q*U[k2][j0]
        nzc=[j for j in range(r2,nk) if H[i][j]]
        if not nzc: continue
        j0=nzc[0]
        if j0!=r2:
            for k2 in range(nr): H[k2][r2],H[k2][j0]=H[k2][j0],H[k2][r2]
            for k2 in range(nk): U[k2][r2],U[k2][j0]=U[k2][j0],U[k2][r2]
        pv.append((i,r2)); r2+=1
    y=[0]*nk
    for i,j in pv:
        s=rhs[i]-sum(H[i][k2]*y[k2] for k2 in range(j))
        if s%H[i][j]: return None
        y[j]=s//H[i][j]
    for i in range(nr):
        if sum(H[i][k2]*y[k2] for k2 in range(nk))!=rhs[i]: return None
    return [sum(U[k2][j]*y[j] for j in range(nk)) for k2 in range(nk)]

def apply_and_score(x, tag, D):
    w=list(v)
    for j,u in enumerate(K): w[u]=x[j]
    av2=L.all_atom_values(w); fe2=L.failing_eqs(av2); s2=L.NEQ-len(fe2)
    print('   %s |D|=%d -> SCORE %d (violate %s)'%(tag,len(D),s2,sorted(EQ[i] for i in D)),flush=True)
    if s2>=39026:
        out='/home/user/integer_solver/solve_lab/agentA_work/A_canon_%d.json'%s2
        json.dump({str(i):str(w[i]) for i in range(L.NVARS)},open(out,'w'))
        print('   SAVED %s'%out,flush=True)
    return s2

ALL=set(range(n))
best=s0
# ---- step 2: targeted drop sets built from the residual atoms' equations ----
resid=[a for a in nz]
cands=set()
poolE=sorted(set(e for a in resid for e in L.atom2eq[a]) & set(EQ))
print('residual-atom equations in region: %d %s'%(len(poolE),poolE),flush=True)
pool=[E2I[e] for e in poolE]
for k in range(1,7):
    for D in itertools.combinations(pool,k): cands.add(frozenset(D))
print('targeted candidate drop-sets: %d'%len(cands),flush=True)
t0=time.time()
for D in sorted(cands,key=len):
    x=int_solve(sorted(ALL-D))
    if x is not None:
        s2=apply_and_score(x,'TARGETED',D)
        if s2>best: best=s2
        if len(D)<=6: break
print('targeted phase done [%.0fs] best=%d'%(time.time()-t0,best),flush=True)
