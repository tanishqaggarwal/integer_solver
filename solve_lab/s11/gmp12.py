"""What is the obstruction, exactly?  Print the certificate coefficients."""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
from gmp7 import allchk, failing, resp_at
from gmp9 import cone_free
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
bd=allchk(base); F=failing(bd)
rows=set(F); knobs=set()
for g in range(10):
    nk=cone_free(rows,3)-knobs
    if not nk and knobs: break
    knobs|=nk
    R2=set(rows)
    for u in sorted(knobs): R2 |= set(resp_at(base,bd,u))
    if R2==rows: break
    rows=R2
rowl=sorted(rows); knobs=sorted(knobs)
idx={a:i for i,a in enumerate(rowl)}
M=[[0]*len(knobs) for _ in rowl]
for j,u in enumerate(knobs):
    for a,x in resp_at(base,bd,u,rowl).items(): M[idx[a]][j]=x
rhs=[(-bd[a])%P for a in rowl]
m=len(rowl); n=len(knobs)
A=[M[i]+[1 if j==i else 0 for j in range(m)]+[rhs[i]] for i in range(m)]
r=0
for c in range(n):
    pr=next((i for i in range(r,m) if A[i][c]),None)
    if pr is None: continue
    A[r],A[pr]=A[pr],A[r]
    inv=pow(A[r][c],-1,P); A[r]=[x*inv%P for x in A[r]]
    for i in range(m):
        if i!=r and A[i][c]:
            f=A[i][c]; A[i]=[(A[i][k]-f*A[r][k])%P for k in range(n+m+1)]
    r+=1
    if r==m: break
certs=[i for i in range(r,m) if not any(A[i][:n]) and A[i][n+m]]
print(f"{m} rows x {n} knobs, rank {r}, {len(certs)} certificates")
best=min(certs, key=lambda i: sum(1 for j in range(m) if A[i][n+j]))
sup=[(rowl[j], A[best][n+j]) for j in range(m) if A[best][n+j]]
print("smallest certificate:")
for a,c in sup:
    print(f"   a{a}  coeff={c}   residue={'FAIL' if bd[a] else '0'}  eqs={len(L.atom2eq.get(a,{}))}")
# normalise by the first coefficient
c0=sup[0][1]; inv=pow(c0,-1,P)
print("  normalised:", [(a, c*inv%P if (c*inv%P)<10**12 else -( (-c*inv)%P ) if ((-c*inv)%P)<10**12 else str(c*inv%P)[:14]+'..') for a,c in sup])
print("  y . rhs =", (sum(c*rhs[idx[a]] for a,c in sup))%P)
