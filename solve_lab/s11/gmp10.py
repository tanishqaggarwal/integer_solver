"""Diagnose the mod-p closure: which rows obstruct, and can deeper knobs reach them?"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from gmp1 import evalp, forwardp
from gmp7 import allchk, failing, resp_at
from gmp9 import cone_free
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(400000)
FREE=set(u for u in range(L.NVARS) if u not in L.definer)
base=[int(x) for x in json.load(open(os.path.join(HERE,'data','gmp1_state.json')))]
depth=int(sys.argv[1]) if len(sys.argv)>1 else 3
bd=allchk(base); F=failing(bd)
rows=set(F); knobs=set()
for g in range(10):
    nk=cone_free(rows,depth)-knobs
    if not nk and knobs: break
    knobs|=nk
    R2=set(rows)
    for u in sorted(knobs): R2 |= set(resp_at(base,bd,u))
    if R2==rows: break
    rows=R2
    if len(rows)>1500 or len(knobs)>1500: break
rowl=sorted(rows); knobs=sorted(knobs)
idx={a:i for i,a in enumerate(rowl)}
M=[[0]*len(knobs) for _ in rowl]
for j,u in enumerate(knobs):
    for a,x in resp_at(base,bd,u,rowl).items(): M[idx[a]][j]=x
rhs=[(-bd[a])%P for a in rowl]
m=len(rowl); n=len(knobs)
print(f"closure depth {depth}: {m} rows x {n} knobs")
# [M | I | rhs] elimination to expose left-kernel obstructions
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
print("  rank mod p =",r)
certs=[i for i in range(r,m) if not any(A[i][:n]) and A[i][n+m]]
print("  obstructing certificates:",len(certs))
for i in certs:
    sup=[rowl[j] for j in range(m) if A[i][n+j]]
    print(f"    support {len(sup)}: {sup}   equations {len(set().union(*[set(L.atom2eq.get(a,{})) for a in sup]))}")
# rows nobody can move
dead=[a for a in rowl if not any(M[idx[a]])]
print("  rows no knob can move:",len(dead), dead[:20])
print("  of those, currently FAILING:", [a for a in dead if bd[a]])
