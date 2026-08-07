"""Rigorous exhaustive floor for the canonical basin: no violation set of size <= K is
mod-p consistent.  (Same reduction as modpobs.py: g in span{col_i : i in D}.)"""
import sys, time, itertools
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P
K_MAX=int(sys.argv[1]) if len(sys.argv)>1 else 4
v=L.load('/home/user/integer_solver/solve_lab/s10/mod9118_0.json')
fe=L.failing_eqs(L.all_atom_values(v))
A=set(a for e in fe for a in L.eq_atoms[e][2])
K,R,rows=build(v,A); nk=len(K)
NZ=[(e,c,lin) for e,c,lin,hq in rows if not hq and lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0)%P for j in range(nk)] for e,c,lin in NZ]
B=[(-c)%P for e,c,lin in NZ]
M=[[N[i][j] for j in range(nk)]+[1 if k==i else 0 for k in range(n)] for i in range(n)]
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
g=[sum(Wb[j][i]*B[i] for i in range(n))%P for j in range(w)]
cols=[[Wb[j][i] for j in range(w)] for i in range(n)]
# precompute a reduced form of g against each column for speed
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
t0=time.time()
for k in range(1,K_MAX+1):
    hit=None; cnt=0
    for S in itertools.combinations(range(n),k):
        cnt+=1
        if in_span(S): hit=[EQ[i] for i in S]; break
    print('|D|=%d : %d subsets exhausted, mod-p consistent: %s  [%.0fs]'%(
        k,cnt,hit if hit else 'NONE',time.time()-t0),flush=True)
    if hit: break
else:
    print('RIGOROUS: no violation set of size <= %d is mod-p consistent -> every integer'%K_MAX)
    print('knob vector in this region violates >= %d equations.'%(K_MAX+1),flush=True)
