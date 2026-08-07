"""Find the mod-p dual certificate of the copy network's obstruction.

The 386-row response system is rationally consistent but insoluble mod p, so there is a left
vector y with  y^T M == 0 (mod p)  and  y . rhs != 0 (mod p).  Its SUPPORT is exactly the set
of rows that must be involved in any repair -- and any drop set has to intersect it.
"""
import sys, os, json, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
import resp as R
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
src=sys.argv[1]; hops=int(sys.argv[2]) if len(sys.argv)>2 else 2
v=load_raw(src)
BR=[a for a in range(L.NA) if R.av(v,a)!=0]
C=R.candidates(v,BR,hops)
cols={}
for u in C:
    d1,_,_=R.response(v,u,1)
    if not d1: continue
    d2,_,_=R.response(v,u,2)
    if all(d2.get(a,0)==2*d1.get(a,0) for a in set(d1)|set(d2)): cols[u]=d1
ROWS=sorted(set(BR)&R.ISCHK | set().union(*[set(d) for d in cols.values()]))
used=sorted(cols)
M=[[cols[u].get(a,0)%P for u in used] for a in ROWS]
rhs=[(-R.av(v,a))%P for a in ROWS]
m=len(ROWS); n=len(used)
# row-reduce [M | I | rhs] over GF(p) to expose left-kernel rows
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
print(f"system {m}x{n}  rank mod p = {r}")
certs=[i for i in range(r,m) if not any(A[i][:n]) and A[i][n+m]]
print(f"obstructing left-kernel rows: {len(certs)}")
best=None
for i in certs:
    sup=[ROWS[j] for j in range(m) if A[i][n+j]]
    cost=len(set().union(*[set(L.atom2eq.get(a,{})) for a in sup]))
    if best is None or len(sup)<len(best[0]): best=(sup,cost)
    if len(sup)<=6:
        print(f"  certificate support {len(sup)}: {sup}  (equations {cost})")
if best:
    sup,cost=best
    print(f"SMALLEST certificate support = {len(sup)} rows, spanning {cost} equations")
    eqc=sorted((len(L.atom2eq.get(a,{})),a) for a in sup)
    print("  rows by #equations:",eqc[:25])
