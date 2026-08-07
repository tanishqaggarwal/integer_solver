"""EXHAUSTIVE minimum-support enumeration of the region code (fixed information set),
then exact HNF integrality of every support of size <= WMAX."""
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
NZ=[(e,c,lin) for e,c,lin in aff if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0) for j in range(nk)] for e,c,lin in NZ]
B=[-c for e,c,lin in NZ]
d0=[MK.v0[u] for u in K]
cur=[EQ[i] for i in range(n) if sum(N[i][j]*d0[j] for j in range(nk))!=B[i]]
print('rows=%d knobs=%d currently violated=%d %s'%(n,nk,len(cur),cur),flush=True)
Nq=[[x%Q for x in r] for r in N]
# --- systematic form with information set I = first rank-nk independent rows
I=[]; M=[]
for i in range(n):
    cand=M+[Nq[i][:]]
    # gaussian rank check
    T=[r[:] for r in cand]; rr=0
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
    if rr==len(cand): I.append(i); M.append(Nq[i][:])
    if len(I)==nk: break
assert len(I)==nk, len(I)
OUT=[i for i in range(n) if i not in I]
print('information set size %d, OUT size %d'%(len(I),len(OUT)),flush=True)
# invert M: find Uc so that codeword for e_t is c_i = sum_j Nq[i][j]*u_j with u=Uc[:,t]
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
# sanity: CW[t] restricted to I is e_t
for t in range(nk):
    assert all(CW[t][I[k]]==(1 if k==t else 0) for k in range(nk))
Aout=[[CW[t][k] for k in OUT] for t in range(nk)]   # nk x |OUT|
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
sup2w={}
t0=time.time()
for j in range(1,WMAX+1):
    need = nout-(WMAX-j)          # number of OUT coords forced to zero
    if need>nout: continue
    cnt=0
    for J in itertools.combinations(range(nk), j):
        AJ=[[Aout[t][k] for k in range(nout)] for t in J]
        for S in itertools.combinations(range(nout), max(need,0)):
            sub=[[AJ[a][s] for a in range(j)] for s in S]      # |S| x j  (rows=S, cols=J)
            ks=kernel(sub, j)
            for u in ks:
                if all(x==0 for x in u): continue
                cw=[0]*n
                for a,t in enumerate(J):
                    if u[a]:
                        for i in range(n): cw[i]=(cw[i]+u[a]*CW[t][i])%Q
                sup=frozenset(EQ[i] for i in range(n) if cw[i])
                if 0<len(sup)<=WMAX: sup2w[sup]=min(sup2w.get(sup,99),len(sup))
        cnt+=1
    print('  j=%d done (%d J-sets) supports so far=%d [%.0fs]'%(j,cnt,len(sup2w),time.time()-t0),flush=True)
print('EXHAUSTIVE: %d supports of size <= %d'%(len(sup2w),WMAX),flush=True)
print('   size histogram:',sorted(collections.Counter(len(s) for s in sup2w).items()),flush=True)
json.dump([sorted(s) for s in sup2w],open('/home/user/integer_solver/solve_lab/agentA_work/supports.json','w'))
