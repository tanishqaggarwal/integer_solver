"""Extract the conserved-obstruction certificate: vectors c over the 27 residual eqs with
c.J = 0 (mod p) for ALL free inputs but c.r != 0. These fixed mod-p values are the wall."""
import heal_harness as H
from jac_lib import D, freeidx
p=H.p
FAILS=[2554,6816,8124,8680,9421,12231,12270,12350,14584,22044,29125]
RIP=[697,1985,5225,10815,16048,17784,17801,22402,23667,24721,27124,28737,29638,29959,35935,37431]
EQS=FAILS+RIP
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
vd=[None]*H.NVARS
for j in H.freeinp: vd[j]=D(H.val[j],{freeidx[j]:1})
ns={'v':vd,'__builtins__':{}}
for k,t in enumerate(H.order):
    r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
rows=[]; rvec=[]
for i in EQS:
    rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
    if isinstance(rr,D): rows.append(dict(rr.g)); rvec.append(rr.v)
    else: rows.append({}); rvec.append(rr%p)
nr=len(rows)
cols=sorted(set(c for g in rows for c in g))
cidx={c:k for k,c in enumerate(cols)}
nc=len(cols)
# Build J (nr x nc) and find LEFT null space: c such that c^T J = 0.  => J^T c = 0, c in null(J^T)
# Work over GF(p). Represent J^T as nc x nr. Solve for c (length nr).
import numpy as np
# dense mod-p elimination on J^T (nc rows, nr cols) to get null space basis of dimension nr-rank
def inv(a): return pow(a%p,p-2,p)
# Build matrix Jt: rows indexed by cols(nc), columns indexed by eqs(nr)
Jt=[[0]*nr for _ in range(nc)]
for ei,g in enumerate(rows):
    for c,v in g.items(): Jt[cidx[c]][ei]=v%p
# Row-reduce Jt (nc x nr); null space (over c-space of dim nr) = solutions to Jt c = 0
A=[row[:] for row in Jt]
m=nc; n=nr
pivcol=[]; r=0
for col in range(n):
    pr=-1
    for rr in range(r,m):
        if A[rr][col]%p!=0: pr=rr; break
    if pr<0: continue
    A[r],A[pr]=A[pr],A[r]
    iv=inv(A[r][col]); A[r]=[(x*iv)%p for x in A[r]]
    for rr in range(m):
        if rr!=r and A[rr][col]%p!=0:
            f=A[rr][col]; A[rr]=[(A[rr][k]-f*A[r][k])%p for k in range(n)]
    pivcol.append(col); r+=1
    if r==m: break
rank=r; free=[c for c in range(n) if c not in pivcol]
print(f'27 residual eqs; J has rank {rank}; left-null (conserved) dim = {n-rank}')
# build null basis
basis=[]
for fc in free:
    c=[0]*n; c[fc]=1
    for ri,pc in enumerate(pivcol): c[pc]=(-A[ri][fc])%p
    basis.append(c)
print(f'extracted {len(basis)} conserved functionals')
for bi,c in enumerate(basis):
    # verify c^T J = 0
    ok=all(sum(c[ei]*rows[ei].get(col,0) for ei in range(nr))%p==0 for col in cols)
    val=sum(c[ei]*rvec[ei] for ei in range(nr))%p
    supp=[EQS[ei] for ei in range(nr) if c[ei]%p!=0]
    print(f'  functional #{bi}: c.J==0 mod p? {ok};  c.residual = {val}  ({"NONZERO -> obstruction" if val!=0 else "zero"})')
    print(f'      support eqs: {supp}')
