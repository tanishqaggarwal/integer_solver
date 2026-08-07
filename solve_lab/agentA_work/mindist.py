"""Minimum support weight of the region code, the RIGOROUS way.
C = column space of N (length n, dim nk).  H = a basis of the left kernel of N
(dim n-nk).  Then  min support of C = min number of LINEARLY DEPENDENT COLUMNS of H.
Exhaustive for sizes 1..4 (and 5 if it fits); plus the direct knob-column supports."""
import sys, json, time, itertools, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from regsolve2 import build
P=env.P; Q=(1<<61)-1
path=sys.argv[1] if len(sys.argv)>1 else '/home/user/integer_solver/solve_lab/s10/mod9118_0.json'
MAXK=int(sys.argv[2]) if len(sys.argv)>2 else 4
v=L.load(path); av=L.all_atom_values(v); fe=L.failing_eqs(av); s0=L.NEQ-len(fe)
A=set(a for e in fe for a in L.eq_atoms[e][2])
K,R,rows=build(v,A); nk=len(K)
good=[(e,c,lin) for e,c,lin,hq in rows if not hq]
NZ=[(e,c,lin) for e,c,lin in good if lin]
n=len(NZ); EQ=[e for e,_,_ in NZ]
N=[[lin.get(j,0)%Q for j in range(nk)] for e,c,lin in NZ]
print('%s: region rows n=%d knobs nk=%d  (dual dim = %d)'%(path.split('/')[-1],n,nk,n-nk),flush=True)
# direct: support of each knob column
colsup=collections.Counter()
for j in range(nk):
    s=sum(1 for i in range(n) if N[i][j])
    colsup[s]+=1
print('knob-column support sizes (upper bounds on d):',sorted(colsup.items())[:10],flush=True)
print('   smallest single-knob support = %d'%min(s for s in colsup.elements()),flush=True)
# left kernel H of N  (vectors lam with lam^T N = 0), dimension n-nk
M=[[N[i][j] for j in range(nk)]+[1 if k==i else 0 for k in range(n)] for i in range(n)]
r=0; piv=[]
for c in range(nk):
    pr=None
    for i in range(r,n):
        if M[i][c]: pr=i;break
    if pr is None: continue
    M[r],M[pr]=M[pr],M[r]
    inv=pow(M[r][c],Q-2,Q); M[r]=[x*inv%Q for x in M[r]]
    for i in range(n):
        if i!=r and M[i][c]:
            f=M[i][c]; M[i]=[(a-f*b)%Q for a,b in zip(M[i],M[r])]
    piv.append(c); r+=1
print('rank(N) mod Q = %d'%r,flush=True)
H=[M[i][nk:] for i in range(r,n)]        # (n-r) x n : rows are left-kernel vectors
m=len(H)
print('parity-check H is %d x %d'%(m,n),flush=True)
cols=[[H[t][i] for t in range(m)] for i in range(n)]
zero=[i for i in range(n) if not any(cols[i])]
print('zero columns of H (weight-1 codewords): %d %s'%(len(zero),[EQ[i] for i in zero][:10]),flush=True)
def rank(vs):
    T=[x[:] for x in vs]; rr=0
    for c in range(m):
        pr=None
        for i in range(rr,len(T)):
            if T[i][c]: pr=i;break
        if pr is None: continue
        T[rr],T[pr]=T[pr],T[rr]
        inv=pow(T[rr][c],Q-2,Q); T[rr]=[x*inv%Q for x in T[rr]]
        for i in range(len(T)):
            if i!=rr and T[i][c]:
                f=T[i][c]; T[i]=[(a-f*b)%Q for a,b in zip(T[i],T[rr])]
        rr+=1
        if rr==len(T): break
    return rr
t0=time.time(); found=None
for k in range(2,MAXK+1):
    cnt=0
    for S in itertools.combinations(range(n),k):
        if rank([cols[i] for i in S])<k:
            found=(k,[EQ[i] for i in S]); break
        cnt+=1
    print('  size %d: %d subsets checked, dependent found: %s  [%.0fs]'%(
        k,cnt,found if found and found[0]==k else 'none',time.time()-t0),flush=True)
    if found: break
if not found:
    print('EXHAUSTIVE: no %d or fewer dependent columns -> minimum support weight >= %d'%(MAXK,MAXK+1),flush=True)
