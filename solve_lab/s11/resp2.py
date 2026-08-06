"""Diagnose the ripple-response system: rational rank, mod-p solvability, obstructing rows,
and what the non-affine columns actually are."""
import sys, os, json, time, collections
from fractions import Fraction
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
from zsolve import solve_int
from ip14 import gf_solve
import resp as R
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)

src=sys.argv[1]; hops=int(sys.argv[2])
v=load_raw(src)
BR=[a for a in range(L.NA) if R.av(v,a)!=0]
C=R.candidates(v,BR,hops)
print(f"{os.path.basename(src)} broken={BR}  candidates={len(C)}")
cols={}; nonlin={}
for u in C:
    d1,_,_=R.response(v,u,1)
    if not d1: continue
    d2,_,_=R.response(v,u,2)
    if all(d2.get(a,0)==2*d1.get(a,0) for a in set(d1)|set(d2)): cols[u]=d1
    else: nonlin[u]=(d1,d2)
print(f"affine={len(cols)} nonaffine={len(nonlin)}")
vals=collections.Counter('bool' if v[u] in (0,1) else ('P' if v[u]==P else 'other') for u in nonlin)
print("  non-affine columns by current value:",dict(vals))
ROWS=sorted(set(BR)&R.ISCHK | set().union(*[set(d) for d in cols.values()]))
used=sorted(cols)
M=[[cols[u].get(a,0) for u in used] for a in ROWS]
rhs=[-R.av(v,a) for a in ROWS]
print(f"system {len(ROWS)}x{len(used)}")
t0=time.time()
xg=gf_solve([[c%P for c in r] for r in M],[c%P for c in rhs],P)
print(f"  solvable mod p: {xg is not None} ({time.time()-t0:.0f}s)")
# rational
def rat_solve(M,rhs):
    m=len(M); n=len(M[0])
    A=[[Fraction(M[i][j]) for j in range(n)]+[Fraction(rhs[i])] for i in range(m)]
    r=0; piv=[]
    for c in range(n):
        pr=next((i for i in range(r,m) if A[i][c]),None)
        if pr is None: continue
        A[r],A[pr]=A[pr],A[r]
        inv=1/A[r][c]; A[r]=[x*inv for x in A[r]]
        for i in range(m):
            if i!=r and A[i][c]:
                f=A[i][c]; A[i]=[A[i][k]-f*A[r][k] for k in range(n+1)]
        piv.append(c); r+=1
        if r==m: break
    bad=[i for i in range(r,m) if A[i][n]!=0 and not any(A[i][:n])]
    return r,piv,A,bad
t0=time.time()
rk,piv,A,bad=rat_solve(M,rhs)
print(f"  rational rank={rk} of {min(len(ROWS),len(used))}; inconsistent rows: {len(bad)} ({time.time()-t0:.0f}s)")
if xg is not None and not bad:
    print("  => rationally & mod-p solvable; the obstruction is INTEGRALITY (denominators)")
# which rows can't be reached at all?
reach=set()
for d in cols.values(): reach|=set(d)
print("  broken check rows with NO column touching them:", [a for a in ROWS if R.av(v,a)!=0 and a not in reach])
