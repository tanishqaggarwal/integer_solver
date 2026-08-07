import sys, json, itertools; sys.path.insert(0,'.')
from fractions import Fraction as F
M=json.load(open('/home/user/integer_solver/solve_lab/agentA_work/model22.json'))
K=M['K']; rows=M['rows']
aff=[r for r in rows if not r['quad']]
n=len(aff); m=len(K)
Lin=[[int(x) for x in r['lin']] for r in aff]
C=[int(r['c']) for r in aff]
EQ=[r['eq'] for r in aff]

def rref(mat):
    mat=[row[:] for row in mat]; nr=len(mat); nc=len(mat[0])
    piv=[]; r=0
    for c in range(nc):
        pr=None
        for i in range(r,nr):
            if mat[i][c]!=0: pr=i;break
        if pr is None: continue
        mat[r],mat[pr]=mat[pr],mat[r]
        pv=mat[r][c]
        mat[r]=[F(x,1)/pv for x in mat[r]]
        for i in range(nr):
            if i!=r and mat[i][c]!=0:
                f=mat[i][c]; mat[i]=[a-f*b for a,b in zip(mat[i],mat[r])]
        piv.append(c); r+=1
        if r==nr: break
    return mat,piv

# rank of lin, and of augmented
_,p1=rref(Lin)
_,p2=rref([Lin[i]+[C[i]] for i in range(n)])
print('rows=%d knobs=%d rank(lin)=%d rank(aug)=%d'%(n,m,len(p1),len(p2)))
# left kernel of Lin: lambda with lambda^T Lin = 0
# compute by rref of Lin^T augmented with identity
T=[[F(Lin[j][i]) for j in range(n)]+[F(1) if j==k else F(0) for k in range(n)] for i in range(m)]
Tr,piv=rref(T)
# rows of Tr whose first n entries are zero give left-kernel vectors
ker=[r[n:] for r in Tr if all(x==0 for x in r[:n])]
# also rows beyond rank
print('left-kernel dim (from rref) =', n-len(p1))
# proper computation: nullspace of Lin^T (m x n) -> vectors in R^n
def nullspace(mat):
    R,piv=rref(mat); nc=len(mat[0]); free=[c for c in range(nc) if c not in piv]
    basis=[]
    for fc in free:
        vec=[F(0)]*nc; vec[fc]=F(1)
        for i,pc in enumerate(piv): vec[pc]=-R[i][fc]
        basis.append(vec)
    return basis
LT=[[F(Lin[j][i]) for j in range(n)] for i in range(m)]
KER=nullspace(LT)
print('left kernel dim =',len(KER))
vals=[sum(k[i]*C[i] for i in range(n)) for k in KER]
nz=[i for i,x in enumerate(vals) if x!=0]
print('certificates (lambda with lambda.c != 0): %d of %d basis vectors'%(len(nz),len(KER)))
# reduce: the functional f on the kernel has rank 0 or 1
print('rank of f on kernel =', 0 if not nz else 1)
# find the support structure: pick basis so that only one has f!=0
if nz:
    i0=nz[0]; base=KER[i0]; v0=vals[i0]
    newk=[]
    for i,k in enumerate(KER):
        if i==i0: continue
        f=vals[i]
        newk.append([a-F(f,1)/v0*b for a,b in zip(k,base)])
    print('certificate vector support (rows that must lose >=1):')
    sup=[EQ[i] for i in range(n) if base[i]!=0]
    print('  ',sup)
    print('   coeffs:', {EQ[i]:str(base[i]) for i in range(n) if base[i]!=0})
