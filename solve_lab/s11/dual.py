"""DUAL CERTIFICATE of the p-obstruction.

The primal says: M x = rhs is unsolvable over Z, obstruction = one factor of p.
The dual says: there is y with  y^T M = 0 (mod p)  but  y^T rhs != 0 (mod p).
Any feasible subsystem must DROP at least one row in supp(y).  So the sparsest such y
is a lower bound on -- and a blueprint for -- the cheapest sacrifice.
"""
import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
from ip8 import build
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); sys.set_int_max_str_digits(400000)

def left_kernel_modp(M, q):
    """basis of {y : y^T M = 0 mod q}, via row reduction of M^T augmented with identity"""
    m=len(M); n=len(M[0])
    # work with rows of M; augment each row with e_i
    A=[[M[i][j]%q for j in range(n)]+[1 if k==i else 0 for k in range(m)] for i in range(m)]
    r=0
    for c in range(n):
        pr=None
        for i in range(r,m):
            if A[i][c]: pr=i; break
        if pr is None: continue
        A[r],A[pr]=A[pr],A[r]
        inv=pow(A[r][c],-1,q)
        A[r]=[x*inv%q for x in A[r]]
        for i in range(m):
            if i!=r and A[i][c]:
                f=A[i][c]
                A[i]=[(A[i][k]-f*A[r][k])%q for k in range(n+m)]
        r+=1
        if r==m: break
    ker=[A[i][n:] for i in range(r,m)]
    return ker, r

LAB=os.path.join(HERE,'..')
src=sys.argv[1] if len(sys.argv)>1 else os.path.join(LAB,'best','new_instance_partial_39026.json')
v=load_raw(src)
print("===", os.path.basename(src))
v,FAIL,used,M,rhs,nf=build(v)
m=len(M)
ker,rank=left_kernel_modp(M,P)
print(f"  rank(M) mod p = {rank}; left-kernel dim = {len(ker)} (of {m} rows)")
# which kernel vectors detect the obstruction?
det=[]
for y in ker:
    s=sum(y[i]*rhs[i] for i in range(m))%P
    if s: det.append(y)
print(f"  kernel vectors with y.rhs != 0 mod p : {len(det)} of {len(ker)}")
if det:
    supps=sorted((sum(1 for t in y if t), i) for i,y in enumerate(det))
    print(f"  sparsest detecting certificate support: {supps[0][0]} rows")
    y=det[supps[0][1]]
    sup=[i for i,t in enumerate(y) if t]
    print(f"  its support (row indices): {sup[:40]}")
    print(f"  of those, FAILING rows: {[i for i in sup if i<nf]}  (nf={nf})")
